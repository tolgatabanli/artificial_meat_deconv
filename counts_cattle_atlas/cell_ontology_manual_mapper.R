## ------------------------------------------------------------------
## Interactive Cell Ontology mapping tool -- standalone script
## Everything is read from scratch: no pre-existing objects required.
## Edit the CONFIG block below, then source/run this file.
## ------------------------------------------------------------------

library(shiny)
library(DT)
library(stringdist)
library(dplyr)
library(ontologyIndex)
library(Seurat)

## ------------------------------------------------------------------
## 0. CONFIG -- edit these paths
## ------------------------------------------------------------------
setwd("/nfs/home/students/t.tabanli/artifical_meat_deconv/counts_cattle_atlas")
CATTLE_RDS <- "sampled_cattle_atlas.RDS"     # path to the Seurat object (.rds)
OBO_PATH   <- "cl-basic.obo"         # local cache path for the ontology
OUT_CSV    <- "celltype_to_CL_map.csv"
N_CANDIDATES <- 5

## ------------------------------------------------------------------
## 1a. Download the Cell Ontology if not already cached locally
## ------------------------------------------------------------------
if (!file.exists(OBO_PATH)) {
  download.file(
    "http://purl.obolibrary.org/obo/cl/cl-basic.obo",
    OBO_PATH, mode = "wb"
  )
}

cl <- get_ontology(
  OBO_PATH,
  propagate_relationships = c("is_a", "part_of"),
  extract_tags = "everything"
)
active_ids <- names(cl$name)[!cl$obsolete]

clean_syn <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(character(0))
  gsub('^"(.*?)".*$', "\\1", x)
}

cl_df <- data.frame(cl_id = active_ids, term = cl$name[active_ids],
                    source = "name", stringsAsFactors = FALSE)
syn_df <- lapply(active_ids, function(id) {
  syns <- clean_syn(cl$synonym[[id]])
  if (length(syns) == 0) return(NULL)
  data.frame(cl_id = id, term = syns, source = "synonym",
             stringsAsFactors = FALSE)
}) %>% bind_rows()

cl_df_clean <- bind_rows(cl_df, syn_df) %>%
  mutate(term_clean = tolower(trimws(term)))

## ------------------------------------------------------------------
## 1b. Load the cattle atlas and pull tissue / cell_type metadata
## ------------------------------------------------------------------
cattle <- readRDS(CATTLE_RDS)
cattle$cell_type <- trimws(cattle$cell_type)
cattle$cell_type <- gsub("'", "", cattle$cell_type) # Henle's loop causes issues in column parse
cattle$cell_type <- gsub("/", "_", cattle$cell_type)

df_cattle <- data.frame(
  tissue   = cattle@meta.data$tissue,
  original = cattle@meta.data$cell_type,
  stringsAsFactors = FALSE
) %>%
  distinct()

rm(cattle)   # metadata extracted; drop the (likely large) object from memory

## ------------------------------------------------------------------
## 2. Build review units: one row per unique cell type, tissues listed
##    for context (adjust to one row per tissue+celltype pair if you
##    genuinely need tissue-specific mappings)
## ------------------------------------------------------------------
units <- df_cattle %>%
  group_by(original) %>%
  summarise(tissues = paste(sort(unique(tissue)), collapse = ", "),
            .groups = "drop") %>%
  arrange(original)

## ------------------------------------------------------------------
## 3. Precompute top-N candidates per unit
## ------------------------------------------------------------------
precompute_candidates <- function(ct, n = N_CANDIDATES) {
  dists <- stringdist(tolower(ct), tolower(cl_df_clean$term_clean), method = "jw")
  ord <- order(dists)[1:n]
  data.frame(
    cl_id   = cl_df_clean$cl_id[ord],
    term    = cl_df_clean$term_clean[ord],
    dist    = round(dists[ord], 3),
    stringsAsFactors = FALSE
  )
}

## ------------------------------------------------------------------
## 4. Resume support: load existing map if present
## ------------------------------------------------------------------
if (file.exists(OUT_CSV)) {
  existing <- read.csv(OUT_CSV, stringsAsFactors = FALSE)
} else {
  existing <- data.frame(
    original = character(), tissues = character(),
    cl_id = character(), cl_name = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
}

todo <- units %>% filter(!original %in% existing$original)

## ------------------------------------------------------------------
## 5. Shiny app
## ------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Cell type -> Cell Ontology mapping"),
  fluidRow(
    column(12,
           div(style = "background:#f0f0f0; padding:10px; margin-bottom:10px;",
               textOutput("progress"))
    )
  ),
  fluidRow(
    column(6,
           h4("Current cell type"),
           verbatimTextOutput("current_label"),
           h5("Top candidates (click to select)"),
           uiOutput("candidate_buttons"),
           hr(),
           actionButton("skip_btn", "Skip / unmapped", class = "btn-warning"),
           actionButton("manual_btn", "Enter CL ID manually", class = "btn-secondary")
    ),
    column(6,
           h4("Search full ontology"),
           textInput("search_box", "Search term_clean (regex ok):", ""),
           DTOutput("search_table")
    )
  ),
  hr(),
  fluidRow(
    column(12,
           actionButton("prev_btn", "<< Previous"),
           actionButton("save_quit", "Save & Quit", class = "btn-danger"),
           textOutput("selected_preview")
    )
  )
)

server <- function(input, output, session) {
  idx <- reactiveVal(1)
  results <- reactiveVal(existing)
  pending_selection <- reactiveVal(NULL)  # list(cl_id, cl_name, method)
  
  current_unit <- reactive({
    req(idx() <= nrow(todo))
    todo[idx(), ]
  })
  
  candidates <- reactive({
    req(idx() <= nrow(todo))
    precompute_candidates(current_unit()$original)
  })
  
  output$progress <- renderText({
    sprintf("Reviewed: %d already saved | Remaining this session: %d / %d",
            nrow(existing), nrow(todo) - idx() + 1, nrow(todo))
  })
  
  output$current_label <- renderText({
    if (idx() > nrow(todo)) return("All done -- click Save & Quit.")
    sprintf("%s\n(tissues: %s)", current_unit()$original, current_unit()$tissues)
  })
  
  output$candidate_buttons <- renderUI({
    req(idx() <= nrow(todo))
    cands <- candidates()
    lapply(seq_len(nrow(cands)), function(i) {
      actionButton(
        inputId = paste0("cand_", i),
        label = sprintf("[%.3f] %s (%s)", cands$dist[i], cands$term[i], cands$cl_id[i]),
        style = "display:block; margin-bottom:4px; width:100%; text-align:left;"
      )
    })
  })
  
  observe({
    lapply(seq_len(N_CANDIDATES), function(i) {
      observeEvent(input[[paste0("cand_", i)]], {
        cands <- candidates()
        if (i <= nrow(cands)) {
          pending_selection(list(cl_id = cands$cl_id[i], cl_name = cands$term[i],
                                 method = "top5"))
          commit_current()
        }
      }, ignoreInit = TRUE)
    })
  })
  
  output$search_table <- renderDT({
    req(input$search_box != "")
    hits <- cl_df_clean[grepl(input$search_box, cl_df_clean$term_clean,
                              ignore.case = TRUE), ]
    datatable(hits, selection = "single", options = list(pageLength = 8))
  })
  
  observeEvent(input$search_table_rows_selected, {
    req(input$search_box != "")
    hits <- cl_df_clean[grepl(input$search_box, cl_df_clean$term_clean,
                              ignore.case = TRUE), ]
    sel <- input$search_table_rows_selected
    if (length(sel) == 1) {
      row <- hits[sel, ]
      pending_selection(list(cl_id = row$cl_id, cl_name = row$term_clean,
                             method = "search"))
      commit_current()
    }
  })
  
  observeEvent(input$manual_btn, {
    showModal(modalDialog(
      title = "Enter CL ID manually",
      textInput("manual_id", "CL ID (e.g. CL:0000542):"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("manual_confirm", "Confirm")
      )
    ))
  })
  
  observeEvent(input$manual_confirm, {
    id <- trimws(input$manual_id)
    nm <- cl_df_clean$term_clean[cl_df_clean$cl_id == id]
    pending_selection(list(cl_id = id, cl_name = ifelse(length(nm), nm[1], NA),
                           method = "manual"))
    removeModal()
    commit_current()
  })
  
  observeEvent(input$skip_btn, {
    pending_selection(list(cl_id = NA_character_, cl_name = NA_character_,
                           method = "unmapped"))
    commit_current()
  })
  
  commit_current <- function() {
    sel <- pending_selection()
    req(sel)
    u <- current_unit()
    new_row <- data.frame(
      original = u$original, tissues = u$tissues,
      cl_id = sel$cl_id, cl_name = sel$cl_name,
      notes = "",
      stringsAsFactors = FALSE
    )
    updated <- results() %>% filter(original != u$original) %>% bind_rows(new_row)
    results(updated)
    write.csv(updated, OUT_CSV, row.names = FALSE)   # save immediately, every step
    pending_selection(NULL)
    idx(idx() + 1)
  }
  
  observeEvent(input$prev_btn, {
    if (idx() > 1) idx(idx() - 1)
  })
  
  observeEvent(input$save_quit, {
    write.csv(results(), OUT_CSV, row.names = FALSE)
    stopApp(invisible(results()))
  })
  
  output$selected_preview <- renderText({
    sprintf("Saved so far: %d rows -> %s", nrow(results()), OUT_CSV)
  })
}

shinyApp(ui, server)