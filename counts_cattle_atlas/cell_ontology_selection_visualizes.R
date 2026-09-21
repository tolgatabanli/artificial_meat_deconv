## ------------------------------------------------------------------
## Ontology graph explorer -- interactive selection of annotation classes
##
## Click a node to TOGGLE it as a final annotation class. Selecting a node
## collapses its whole subtree into it: descendants are hidden and every
## mapped cell type below it is assigned to it. Click again to expand.
## The set of selected nodes IS your new annotation vocabulary.
##
## The network is rendered ONCE and mutated via visNetworkProxy, so zoom
## and pan are never lost when toggling.
## ------------------------------------------------------------------

library(ontologyIndex)
library(igraph)
library(visNetwork)
library(dplyr)
library(DT)
library(shiny)

OBO_PATH <- "cl-basic.obo"
MAP_CSV  <- "celltype_to_CL_map.csv"   # output of the mapping app
SEL_FILE <- "selected_classes.txt"     # selection is saved/restored here
ROOT_ID  <- "CL:0000000"               # "cell", the CL root

LEVEL_SEP <- 300    # horizontal gap between ontology levels (LR layout)
NODE_SEP  <- 75     # vertical gap between siblings

## ------------------------------------------------------------------
## 1. Load ontology + curated mapping
## ------------------------------------------------------------------
if (!file.exists(OBO_PATH)) {
  download.file("http://purl.obolibrary.org/obo/cl/cl-basic.obo", OBO_PATH, mode = "wb")
}
cl <- get_ontology(OBO_PATH, propagate_relationships = c("is_a", "part_of"),
                   extract_tags = "everything")

map <- read.csv(MAP_CSV, stringsAsFactors = FALSE) %>%
  filter(!is.na(cl_id), cl_id != "")
stopifnot(all(map$cl_id %in% names(cl$name)))   # catch typo'd manual CL IDs early

leaf_ids <- unique(map$cl_id)

## ------------------------------------------------------------------
## 2. Induced subgraph over leaves + their ancestors
## ------------------------------------------------------------------
all_ids <- unique(c(ROOT_ID, unlist(lapply(leaf_ids, function(id) get_ancestors(cl, id)))))

edges_df <- do.call(rbind, lapply(all_ids, function(id) {
  parents <- intersect(cl$parents[[id]], all_ids)
  if (length(parents) == 0) return(NULL)
  data.frame(parent = parents, child = id, stringsAsFactors = FALSE)
}))

g <- graph_from_data_frame(edges_df, directed = TRUE,
                           vertices = data.frame(name = all_ids, stringsAsFactors = FALSE))

depth_of <- distances(g, v = ROOT_ID, to = V(g), mode = "out")[1, ]
names(depth_of) <- V(g)$name

desc_list <- lapply(all_ids, function(id) setdiff(names(subcomponent(g, id, mode = "out")), id))
names(desc_list) <- all_ids

anc_list <- lapply(all_ids, function(id) intersect(get_ancestors(cl, id), all_ids))
names(anc_list) <- all_ids

children_of <- lapply(all_ids, function(id) intersect(cl$children[[id]], all_ids))
names(children_of) <- all_ids

n_below <- vapply(all_ids, function(id) {
  length(intersect(leaf_ids, c(id, desc_list[[id]])))
}, integer(1))
names(n_below) <- all_ids

## all_ids is ancestor-closed, so subgraph parents == full CL parents here
n_parents <- vapply(all_ids, function(id) length(intersect(cl$parents[[id]], all_ids)),
                    integer(1))
names(n_parents) <- all_ids

cat(sprintf("Subgraph: %d nodes, %d of them mapped cell types, max depth %d\n",
            length(all_ids), length(leaf_ids), max(depth_of, na.rm = TRUE)))

## ------------------------------------------------------------------
## 3. Helpers
## ------------------------------------------------------------------

wrap_label <- function(x, width = 20) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

visible_nodes <- function(selected, compress_chains = TRUE) {
  hidden <- unique(unlist(desc_list[selected]))
  ## a selected node is ALWAYS shown, even if it sits under another selected
  ## node -- otherwise a nested class exists in the export but not on screen
  vis <- union(setdiff(all_ids, hidden), selected)
  if (compress_chains) {
    repeat {
      drop <- vapply(vis, function(v) {
        if (v == ROOT_ID || v %in% selected || v %in% leaf_ids) return(FALSE)
        length(intersect(children_of[[v]], vis)) == 1
      }, logical(1))
      if (!any(drop)) break
      vis <- vis[!drop]
    }
  }
  vis
}

nearest_visible_parents <- function(v, vis) {
  out <- character(0); seen <- character(0)
  queue <- intersect(cl$parents[[v]], all_ids)
  while (length(queue)) {
    p <- queue[1]; queue <- queue[-1]
    if (p %in% seen) next
    seen <- c(seen, p)
    if (p %in% vis) out <- c(out, p)
    else queue <- c(queue, intersect(cl$parents[[p]], all_ids))
  }
  unique(out)
}

display_edges <- function(vis) {
  ed <- do.call(rbind, lapply(vis, function(v) {
    ps <- nearest_visible_parents(v, vis)
    if (!length(ps)) return(NULL)
    data.frame(from = ps, to = v, stringsAsFactors = FALSE)
  }))
  if (is.null(ed)) ed <- data.frame(from = character(0), to = character(0),
                                    stringsAsFactors = FALSE)
  ed$id <- paste0(ed$from, "->", ed$to)   # stable ids so proxy can add/remove
  ed
}

## Layout computed in R: layer from ontology depth (compacted to consecutive
## ranks), within-layer ordering from Sugiyama crossing minimisation.
compute_positions <- function(vis, ed) {
  layer <- as.integer(factor(depth_of[vis]))          # 1..k, gaps removed
  names(layer) <- vis
  
  sg <- graph_from_data_frame(
    ed[, c("from", "to"), drop = FALSE], directed = TRUE,
    vertices = data.frame(name = vis, stringsAsFactors = FALSE)
  )
  lay <- layout_with_sugiyama(sg, layers = layer[V(sg)$name])$layout
  
  data.frame(
    id = vis,
    x  = (layer[vis] - 1) * LEVEL_SEP,
    y  = lay[, 1] * NODE_SEP,
    stringsAsFactors = FALSE
  )
}

## Flat scalar columns only -- nested color.*/font.* do not survive
## visUpdateNodes reliably, so colour and font go in as plain strings.
build_nodes <- function(selected, vis, pos) {
  pos_map_x <- setNames(pos$x, pos$id)
  pos_map_y <- setNames(pos$y, pos$id)
  
  ## selected terms that sit underneath another selected term
  shadowed <- selected[vapply(selected, function(s) {
    length(intersect(anc_list[[s]], setdiff(selected, s))) > 0
  }, logical(1))]
  
  data.frame(id = all_ids, stringsAsFactors = FALSE) %>%
    mutate(
      is_vis    = id %in% vis,
      is_sel    = id %in% selected,
      is_shad   = id %in% shadowed,
      is_leaf   = id %in% leaf_ids,
      multi_par = n_parents[id] > 1,
      label     = paste0(ifelse(is_shad, "\u26a0 ", ""),
                         ifelse(multi_par, "* ", ""),
                         wrap_label(cl$name[id])),
      shape     = ifelse(is_leaf | is_sel, "box", "ellipse"),
      color     = case_when(
        is_shad            ~ "#6a040f",   # nested selection -- almost certainly unintended
        is_sel & multi_par ~ "#7209b7",   # selected AND multi-parent -- check for ambiguity
        is_sel             ~ "#e63946",
        is_leaf            ~ "#f4a261",
        TRUE               ~ "#a8dadc"
      ),
      font        = ifelse(is_sel, "18px arial #ffffff", "14px arial #000000"),
      borderWidth = ifelse(is_sel, 3, 1),
      hidden      = !is_vis,
      x           = ifelse(is_vis, pos_map_x[id], 0),
      y           = ifelse(is_vis, pos_map_y[id], 0),
      title       = paste0(id, "<br>", cl$name[id],
                           "<br>depth ", depth_of[id],
                           "<br>", n_parents[id], " parent(s)",
                           "<br>", n_below[id], " mapped cell type(s) at/below",
                           ifelse(is_shad,
                                  "<br><b>nested under another selected class</b>", ""),
                           ifelse(multi_par,
                                  "<br><b>multiple CL parents</b>", ""))
    ) %>%
    select(id, label, shape, color, font, borderWidth, hidden, x, y, title)
}

## Assign each mapped cell type to its most specific selected ancestor.
## In a DAG the "deepest" ancestor is not always unique, so the rule is:
## keep the selected ancestors that have no other selected ancestor below
## them; if more than one survives, the assignment is genuinely ambiguous
## and gets flagged rather than silently resolved.
assign_groups <- function(selected) {
  if (!length(selected)) {
    return(map %>% mutate(group_id = NA_character_, group_name = NA_character_,
                          n_candidates = 0L, candidates = NA_character_))
  }
  
  per_leaf <- lapply(leaf_ids, function(lf) {
    hits <- intersect(anc_list[[lf]], selected)
    if (!length(hits)) {
      return(data.frame(cl_id = lf, group_id = NA_character_,
                        n_candidates = 0L, candidates = NA_character_,
                        stringsAsFactors = FALSE))
    }
    ## most specific = no other hit is a descendant of it
    keep <- vapply(hits, function(h) {
      others <- setdiff(hits, h)
      !any(vapply(others, function(h2) h %in% anc_list[[h2]], logical(1)))
    }, logical(1))
    ms <- hits[keep]
    
    data.frame(
      cl_id        = lf,
      group_id     = ms[which.max(depth_of[ms])],   # deterministic tiebreak
      n_candidates = length(ms),
      candidates   = paste(cl$name[ms], collapse = " | "),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  map %>%
    left_join(per_leaf, by = "cl_id") %>%
    mutate(group_name = ifelse(is.na(group_id), NA_character_, cl$name[group_id]))
}

## ------------------------------------------------------------------
## 4. Validating a selection
## ------------------------------------------------------------------

## Nodes in the subgraph reachable by exactly one parent -- selecting only
## from these guarantees a tree-shaped frontier.
single_parent_nodes <- function() all_ids[n_parents[all_ids] == 1]

## Full diagnostic on a candidate selection.
##   $ambiguous_leaves  -- THE blocker: cell types with >1 incomparable class
##   $multi_parent      -- selected nodes with several CL parents (informational)
##   $nested_classes    -- one selected node is an ancestor of another
##   $unassigned        -- cell types with no selected ancestor
check_selection <- function(selected, verbose = TRUE) {
  selected <- intersect(unique(selected), all_ids)
  res <- assign_groups(selected)
  
  ambiguous <- res %>%
    filter(n_candidates > 1) %>%
    select(original, cl_name, n_candidates, candidates)
  
  unassigned <- res %>% filter(n_candidates == 0) %>% select(original, cl_name)
  
  multi_parent <- if (length(selected)) {
    data.frame(
      cl_id   = selected[n_parents[selected] > 1],
      stringsAsFactors = FALSE
    ) %>%
      mutate(cl_name = cl$name[cl_id],
             n_parents = n_parents[cl_id],
             parents = vapply(cl_id, function(i)
               paste(cl$name[intersect(cl$parents[[i]], all_ids)], collapse = " | "),
               character(1)))
  } else data.frame()
  
  nested <- if (length(selected) > 1) {
    do.call(rbind, lapply(selected, function(a) {
      desc_sel <- selected[selected != a & vapply(selected[selected != a],
                                                  function(b) a %in% anc_list[[b]],
                                                  logical(1))]
      if (!length(desc_sel)) return(NULL)
      data.frame(ancestor_class = cl$name[a], descendant_class = cl$name[desc_sel],
                 stringsAsFactors = FALSE)
    }))
  } else NULL
  if (is.null(nested)) nested <- data.frame()
  
  out <- list(
    n_classes        = length(selected),
    ambiguous_leaves = ambiguous,
    multi_parent     = multi_parent,
    nested_classes   = nested,
    unassigned       = unassigned,
    is_clean         = nrow(ambiguous) == 0 && nrow(unassigned) == 0
  )
  
  if (verbose) {
    message(sprintf("Selection: %d classes | %s",
                    out$n_classes,
                    if (out$is_clean) "CLEAN" else "issues found"))
    if (nrow(ambiguous)) {
      message(sprintf("  AMBIGUOUS (%d cell types have competing classes):", nrow(ambiguous)))
      for (i in seq_len(min(10, nrow(ambiguous)))) {
        message(sprintf("    %s -> %s", ambiguous$original[i], ambiguous$candidates[i]))
      }
    }
    if (nrow(unassigned)) {
      message(sprintf("  UNASSIGNED (%d): %s", nrow(unassigned),
                      paste(utils::head(unassigned$original, 8), collapse = " | ")))
    }
    if (nrow(multi_parent)) {
      message(sprintf("  note: %d selected class(es) have multiple CL parents: %s",
                      nrow(multi_parent),
                      paste(multi_parent$cl_name, collapse = " | ")))
    }
    if (nrow(nested)) {
      message(sprintf("  note: %d nested class pair(s) -- descendant wins for its own leaves",
                      nrow(nested)))
    }
  }
  
  invisible(out)
}

## ------------------------------------------------------------------
## 5. Applying a finished selection to real annotations
##
## These work standalone. To load them without launching the app:
##     SKIP_APP <- TRUE; source("celltype_ontology_graph_explorer.R")
## ------------------------------------------------------------------

## ---- antichain checks -------------------------------------------------
## A valid selection is an ANTICHAIN: no selected term is an ancestor of
## another. Nesting is not fatal (assign_groups takes the deepest selected
## ancestor, so the specific class wins) but it is almost always unintended,
## and it can leave a selected class with zero cell types assigned to it.

## every (descendant, ancestor) pair that both sit in the selection
selected_nesting <- function(selected) {
  selected <- unique(intersect(selected, all_ids))
  empty <- data.frame(descendant = character(0), descendant_name = character(0),
                      ancestor = character(0), ancestor_name = character(0),
                      stringsAsFactors = FALSE)
  if (length(selected) < 2) return(empty)
  
  res <- do.call(rbind, lapply(selected, function(s) {
    ans <- intersect(anc_list[[s]], setdiff(selected, s))   # anc_list includes self
    if (!length(ans)) return(NULL)
    data.frame(descendant = s, descendant_name = cl$name[s],
               ancestor = ans, ancestor_name = cl$name[ans],
               stringsAsFactors = FALSE)
  }))
  if (is.null(res)) empty else res
}

is_antichain <- function(selected) nrow(selected_nesting(selected)) == 0

## drop one side of every nested pair
##   keep = "specific" -> discard any selected term that has a selected descendant
##   keep = "general"  -> discard any selected term that has a selected ancestor
prune_to_antichain <- function(selected, keep = c("specific", "general")) {
  keep <- match.arg(keep)
  selected <- unique(intersect(selected, all_ids))
  nest <- selected_nesting(selected)
  if (!nrow(nest)) return(selected)
  drop <- if (keep == "specific") unique(nest$ancestor) else unique(nest$descendant)
  setdiff(selected, drop)
}

## full health report on a selection
validate_selection <- function(selected, verbose = TRUE) {
  selected <- unique(intersect(selected, all_ids))
  nest <- selected_nesting(selected)
  res  <- assign_groups(selected)
  
  assigned_counts <- vapply(selected, function(s) sum(res$group_id == s, na.rm = TRUE),
                            integer(1))
  out <- list(
    n_classes     = length(selected),
    nesting       = nest,
    empty_classes = names(assigned_counts)[assigned_counts == 0],
    unassigned    = res$original[is.na(res$group_id)]
  )
  
  if (verbose) {
    message(sprintf("%d classes selected", out$n_classes))
    if (nrow(nest)) {
      message(sprintf("  WARNING: %d nested pair(s) -- selection is not an antichain:",
                      nrow(nest)))
      for (i in seq_len(min(nrow(nest), 10))) {
        message(sprintf("    '%s' is under '%s'", nest$descendant_name[i], nest$ancestor_name[i]))
      }
      if (nrow(nest) > 10) message(sprintf("    ... (+%d more)", nrow(nest) - 10))
    } else {
      message("  OK: no selected term is an ancestor of another")
    }
    if (length(out$empty_classes)) {
      message("  WARNING: class(es) with 0 cell types assigned: ",
              paste(cl$name[out$empty_classes], collapse = " | "))
    }
    if (length(out$unassigned)) {
      message(sprintf("  %d cell type(s) still unassigned", length(out$unassigned)))
    }
  }
  invisible(out)
}

## Rebuild the leaf -> class table from a bare set of selected CL IDs.
## Use this when you have selected_classes.txt but not the exported CSV.
build_annotation_map <- function(selected_ids, strict = FALSE) {
  selected_ids <- intersect(unique(selected_ids), all_ids)
  nest <- selected_nesting(selected_ids)
  if (nrow(nest)) {
    msg <- sprintf("selection is not an antichain: %d nested pair(s), e.g. '%s' under '%s'",
                   nrow(nest), nest$descendant_name[1], nest$ancestor_name[1])
    if (strict) stop(msg) else warning(msg, " -- deepest selected ancestor wins")
  }
  assign_groups(selected_ids)
}

## Map a vector of existing annotations onto the selected classes.
##   x       : character/factor vector -- original cell type labels OR CL IDs
##   ann_map : data.frame with original / cl_id / group_id / group_name
##             (the exported CSV, or build_annotation_map(...))
apply_annotation_map <- function(x, ann_map,
                                 by         = c("auto", "original", "cl_id"),
                                 value      = c("group_name", "group_id"),
                                 unassigned = NA_character_,
                                 verbose    = TRUE) {
  by    <- match.arg(by)
  value <- match.arg(value)
  
  need <- c("original", "cl_id", "group_id", "group_name")
  absent <- setdiff(need, names(ann_map))
  if (length(absent)) stop("ann_map is missing column(s): ", paste(absent, collapse = ", "))
  
  x <- as.character(x)
  
  if (by == "auto") {
    by <- if (mean(grepl("^CL:[0-9]+$", x), na.rm = TRUE) > 0.5) "cl_id" else "original"
    if (verbose) message("apply_annotation_map(): matching by '", by, "'")
  }
  
  ## collapse the key -> value table, flagging genuine conflicts
  tab <- unique(data.frame(key = as.character(ann_map[[by]]),
                           val = as.character(ann_map[[value]]),
                           stringsAsFactors = FALSE))
  conflicted <- unique(tab$key[duplicated(tab$key)])
  if (length(conflicted)) {
    warning("inconsistent mapping for key(s): ", paste(conflicted, collapse = ", "),
            " -- first occurrence used")
    tab <- tab[!duplicated(tab$key), ]
  }
  
  lut <- setNames(tab$val, tab$key)
  out <- unname(lut[x])
  
  if (verbose) {
    unmatched  <- setdiff(unique(x), names(lut))
    n_no_key   <- sum(x %in% unmatched)
    n_no_class <- sum(!x %in% unmatched & is.na(out))
    message(sprintf("  %d/%d entries assigned | %d with no map entry | %d mapped but unassigned",
                    sum(!is.na(out)), length(x), n_no_key, n_no_class))
    if (length(unmatched)) {
      message("  labels absent from ann_map: ",
              paste(utils::head(unmatched, 10), collapse = " | "),
              if (length(unmatched) > 10) sprintf(" ... (+%d more)", length(unmatched) - 10) else "")
    }
  }
  
  out[is.na(out)] <- unassigned
  out
}

## Convenience wrapper: write the broad annotation straight into a Seurat object
add_broad_annotation <- function(obj, ann_map, from = "cell_type",
                                 to = "cell_type_broad", ...) {
  if (!from %in% colnames(obj@meta.data)) {
    stop("column '", from, "' not found in obj@meta.data")
  }
  obj@meta.data[[to]] <- apply_annotation_map(obj@meta.data[[from]], ann_map, ...)
  obj
}

## Per-class cell counts, with a flag for classes too small to give a stable
## deconvolution signature.
annotation_summary <- function(x, ann_map, min_cells = 100, ...) {
  grp <- apply_annotation_map(x, ann_map, unassigned = "<unassigned>", verbose = FALSE, ...)
  data.frame(class = grp, stringsAsFactors = FALSE) %>%
    count(class, name = "n_cells") %>%
    left_join(
      ann_map %>% filter(!is.na(group_name)) %>%
        count(group_name, name = "n_source_labels"),
      by = c("class" = "group_name")
    ) %>%
    mutate(too_small = n_cells < min_cells) %>%
    arrange(desc(n_cells))
}

## ------------------------------------------------------------------
## 5. Shiny app
## ------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Cell Ontology -- build annotation classes by toggling nodes"),
  fluidRow(
    column(12, div(
      style = "background:#f4f4f4; padding:8px; margin-bottom:8px;",
      strong("Click a node to toggle it as an annotation class."),
      " Selected nodes turn red and absorb everything below them. ",
      textOutput("status", inline = TRUE)
    ))
  ),
  fluidRow(
    column(12, uiOutput("nest_banner")),
    column(12, uiOutput("ambig_banner"))
  ),
  fluidRow(
    column(8,
           checkboxInput("compress", "Hide pass-through nodes (single-child chains)", TRUE),
           actionButton("clear", "Clear selection"),
           actionButton("preset_leaves", "Select all cell types (identity)"),
           actionButton("fit", "Re-fit view"),
           br(), br(),
           visNetworkOutput("net", height = "780px")
    ),
    column(4,
           h4("Selected classes"),
           DTOutput("sel_table"),
           hr(),
           h4("Unassigned cell types"),
           DTOutput("unassigned"),
           hr(),
           downloadButton("dl", "Download annotation map")
    )
  )
)

server <- function(input, output, session) {
  
  ## restore a previous session's selection if one was saved
  init_sel <- character(0)
  if (file.exists(SEL_FILE)) {
    raw <- trimws(readLines(SEL_FILE, warn = FALSE))
    raw <- raw[nzchar(raw)]
    init_sel <- intersect(raw, all_ids)
    if (length(setdiff(raw, all_ids))) {
      warning("dropped ", length(setdiff(raw, all_ids)),
              " saved ID(s) not present in this subgraph")
    }
  }
  
  selected       <- reactiveVal(init_sel)
  rendered_edges <- reactiveVal(character(0))
  
  ## autosave: the selection survives a crash or a closed browser tab
  observeEvent(selected(), {
    writeLines(selected(), SEL_FILE)
  }, ignoreNULL = FALSE)
  
  graph_data <- reactive({
    vis <- visible_nodes(selected(), input$compress)
    ed  <- display_edges(vis)
    pos <- compute_positions(vis, ed)
    list(nodes = build_nodes(selected(), vis, pos), edges = ed)
  })
  
  ## ---- rendered exactly once -------------------------------------
  output$net <- renderVisNetwork({
    gd <- isolate(graph_data())
    rendered_edges(gd$edges$id)
    
    visNetwork(gd$nodes, gd$edges) %>%
      visNodes(widthConstraint = list(minimum = 90, maximum = 150),
               margin = list(top = 8, bottom = 8, left = 10, right = 10)) %>%
      visEdges(arrows = "to", color = list(color = "#b0b0b0"), smooth = FALSE) %>%
      visPhysics(enabled = FALSE) %>%   # positions come from R, nothing drifts
      visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE,
                     navigationButtons = TRUE, tooltipDelay = 150) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
      visEvents(selectNode = "function(nodes){
         Shiny.setInputValue('node_clicked', nodes.nodes[0], {priority: 'event'});
      }")
  })
  
  ## ---- all subsequent changes go through the proxy ----------------
  observeEvent(list(selected(), input$compress), {
    gd  <- graph_data()
    old <- rendered_edges()
    new <- gd$edges$id
    
    proxy <- visNetworkProxy("net")
    
    gone <- setdiff(old, new)
    if (length(gone)) proxy %>% visRemoveEdges(id = gone)
    
    proxy %>% visUpdateNodes(nodes = gd$nodes)
    if (nrow(gd$edges)) proxy %>% visUpdateEdges(edges = gd$edges)
    
    rendered_edges(new)
  }, ignoreInit = TRUE)   # skip the initial render, which already drew everything
  
  observeEvent(input$node_clicked, {
    id  <- input$node_clicked
    cur <- selected()
    if (id %in% cur) {
      selected(setdiff(cur, id))
    } else {
      selected(c(setdiff(cur, desc_list[[id]]), id))
    }
  })
  
  observeEvent(input$clear, selected(character(0)))
  observeEvent(input$preset_leaves, selected(leaf_ids))
  observeEvent(input$fit, visNetworkProxy("net") %>% visFit())
  
  result <- reactive(assign_groups(selected()))
  
  nesting <- reactive(selected_nesting(selected()))
  
  output$nest_banner <- renderUI({
    nest <- nesting()
    if (!nrow(nest)) {
      return(div(style = "background:#d8f3dc; border-left:5px solid #2d6a4f; padding:8px; margin-bottom:8px;",
                 "\u2713 Selection is a valid antichain \u2014 no selected term is an ancestor of another."))
    }
    div(style = "background:#ffe5e5; border-left:5px solid #9d0208; padding:8px; margin-bottom:8px;",
        strong(sprintf("\u26a0 %d nested pair(s): selected terms that are ancestors of other selected terms.",
                       nrow(nest))),
        br(),
        "Cell types below the deeper term go to the deeper term, so the outer class may end up empty. ",
        "Offending nodes are dark red and flagged \u26a0 in the graph.",
        br(), br(),
        actionButton("prune_specific", "Keep the more specific term(s)", class = "btn-danger"),
        actionButton("prune_general",  "Keep the more general term(s)",  class = "btn-warning")
    )
  })
  
  observeEvent(input$prune_specific, {
    selected(prune_to_antichain(selected(), keep = "specific"))
  })
  observeEvent(input$prune_general, {
    selected(prune_to_antichain(selected(), keep = "general"))
  })
  
  ## The real blocker: a cell type reachable from two incomparable classes.
  ## Multi-parent selected nodes are only a *risk factor* for this, so the
  ## check is done per cell type, not per selected node.
  output$ambig_banner <- renderUI({
    res <- result()
    if (!"n_candidates" %in% names(res)) return(NULL)
    amb <- res %>% filter(n_candidates > 1)
    
    mp <- selected()[n_parents[selected()] > 1]
    
    if (!nrow(amb)) {
      msg <- "\u2713 No ambiguous assignments \u2014 every cell type has exactly one most-specific class."
      if (length(mp)) {
        msg <- paste0(msg, sprintf(" (%d selected class(es) have multiple CL parents, marked * \u2014 harmless here.)",
                                   length(mp)))
      }
      return(div(style = "background:#d8f3dc; border-left:5px solid #2d6a4f; padding:8px; margin-bottom:8px;",
                 msg))
    }
    
    div(style = "background:#fff3cd; border-left:5px solid #b45309; padding:8px; margin-bottom:8px;",
        strong(sprintf("\u26a0 %d cell type(s) have competing classes (DAG ambiguity):", nrow(amb))),
        tags$ul(lapply(seq_len(min(8, nrow(amb))), function(i) {
          tags$li(sprintf("%s \u2192 %s", amb$original[i], amb$candidates[i]))
        })),
        if (nrow(amb) > 8) em(sprintf("... and %d more", nrow(amb) - 8)),
        br(),
        "The exported map resolves these by depth, which is arbitrary. ",
        "Deselect one of each competing pair, or pick a common ancestor instead."
    )
  })
  
  output$status <- renderText({
    res <- result()
    n_amb <- if ("n_candidates" %in% names(res)) sum(res$n_candidates > 1, na.rm = TRUE) else 0L
    sprintf("%d classes selected | %d of %d cell types assigned | %d ambiguous | %d single-parent nodes available",
            length(selected()), sum(!is.na(res$group_id)), nrow(res),
            n_amb, length(single_parent_nodes()))
  })
  
  output$sel_table <- renderDT({
    sel <- selected()
    if (!length(sel)) {
      return(datatable(data.frame(note = "nothing selected yet"),
                       rownames = FALSE, options = list(dom = "t")))
    }
    res <- result()
    data.frame(
      class        = cl$name[sel],
      cl_id        = sel,
      n_cell_types = vapply(sel, function(s) sum(res$group_id == s, na.rm = TRUE), integer(1)),
      nested       = ifelse(sel %in% nesting()$descendant | sel %in% nesting()$ancestor,
                            "\u26a0", ""),
      stringsAsFactors = FALSE
    ) %>%
      arrange(n_cell_types == 0, desc(n_cell_types)) %>%
      datatable(rownames = FALSE, options = list(pageLength = 15, dom = "tp"))
  })
  
  output$unassigned <- renderDT({
    result() %>%
      filter(is.na(group_id)) %>%
      select(original, cl_name) %>%
      datatable(rownames = FALSE, options = list(pageLength = 8, dom = "tp"))
  })
  
  output$dl <- downloadHandler(
    filename = function() sprintf("celltype_annotation_classes_%s.csv", Sys.Date()),
    content  = function(file) write.csv(result(), file, row.names = FALSE)
  )
}

if (!isTRUE(get0("SKIP_APP", ifnotfound = FALSE))) {
  shinyApp(ui, server)
}