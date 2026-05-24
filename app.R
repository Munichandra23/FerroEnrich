## ========================== PACKAGES ==========================
suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(DESeq2)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(igraph)
  library(visNetwork)
  library(STRINGdb)
  library(Cairo)
  library(ragg)
  library(bslib)
})

## ---- optional packages flags (used for graceful fallbacks) ---
.has_readxl <- requireNamespace("readxl", quietly = TRUE)
.has_hs     <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
.has_mm     <- requireNamespace("org.Mm.eg.db", quietly = TRUE)
.has_bm     <- requireNamespace("biomaRt", quietly = TRUE)
.has_gp     <- requireNamespace("gprofiler2", quietly = TRUE)
.has_pheat  <- requireNamespace("pheatmap", quietly = TRUE)

## ======================== SMALL HELPERS =======================
`%||%` <- function(a, b) if (!is.null(a)) a else b

strip_ensembl_version <- function(x) {
  x <- as.character(x)
  sub("\\.\\d+$", "", x)
}

detect_species_from_ids <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return("Human")
  id <- toupper(as.character(x))
  if (startsWith(id, "ENSMUSG")) return("Mouse")
  if (startsWith(id, "ENSG"))    return("Human")
  "Human"
}

species_organism_code <- function(sp) {
  if (sp == "Human") "hsapiens" else "mmusculus"
}

## ---- Ensembl -> Symbol mapping (OrgDb -> biomaRt -> g:Profiler) ----
map_ensembl_to_symbol <- function(ensembl_ids, species = c("Human", "Mouse")) {
  species <- match.arg(species)
  
  ids_raw <- as.character(ensembl_ids)
  ids <- strip_ensembl_version(trimws(ids_raw))
  ids[ids == ""] <- NA_character_
  
  try_orgdb <- function(ids, species) {
    if (species == "Human" && .has_hs) {
      orgdb <- org.Hs.eg.db::org.Hs.eg.db
    } else if (species == "Mouse" && .has_mm) {
      orgdb <- org.Mm.eg.db::org.Mm.eg.db
    } else {
      return(rep(NA_character_, length(ids)))
    }
    
    valid <- tryCatch(
      AnnotationDbi::keys(orgdb, keytype = "ENSEMBL"),
      error = function(e) character()
    )
    
    use_ids <- intersect(ids[!is.na(ids)], valid)
    
    mapped <- rep(NA_character_, length(ids))
    names(mapped) <- ids
    
    if (length(use_ids) > 0) {
      sym <- tryCatch(
        AnnotationDbi::mapIds(
          orgdb,
          keys = use_ids,
          keytype = "ENSEMBL",
          column = "SYMBOL",
          multiVals = "first"
        ),
        error = function(e) NULL
      )
      
      if (!is.null(sym) && length(sym) > 0) {
        mapped[names(sym)] <- unname(sym)
      }
    }
    
    unname(mapped)
  }
  
  gene <- try_orgdb(ids, species)
  
  if (FALSE && any(is.na(gene)) && .has_bm) {
    q_ids <- unique(ids[is.na(gene) & !is.na(ids)])
    
    get_mart_safe <- function(species) {
      dataset <- if (species == "Human") "hsapiens_gene_ensembl" else "mmusculus_gene_ensembl"
      mirrors <- c("useast", "uswest", "asia")
      
      for (m in mirrors) {
        mart <- tryCatch(
          biomaRt::useEnsembl(
            biomart = "genes",
            dataset = dataset,
            mirror = m
          ),
          error = function(e) NULL
        )
        if (!is.null(mart)) return(mart)
      }
      
      mart <- tryCatch(
        biomaRt::useEnsembl(
          biomart = "genes",
          dataset = dataset,
          version = 113
        ),
        error = function(e) NULL
      )
      
      mart
    }
    
    if (length(q_ids) > 0) {
      bm <- get_mart_safe(species)
      
      if (!is.null(bm)) {
        sym_attr <- if (species == "Human") "hgnc_symbol" else "mgi_symbol"
        
        q <- tryCatch(
          biomaRt::getBM(
            attributes = c("ensembl_gene_id", sym_attr),
            filters = "ensembl_gene_id",
            values = q_ids,
            mart = bm
          ),
          error = function(e) NULL
        )
        
        if (!is.null(q) && nrow(q) > 0) {
          q <- q[!is.na(q[[sym_attr]]) & q[[sym_attr]] != "", , drop = FALSE]
          q <- q[!duplicated(q$ensembl_gene_id), , drop = FALSE]
          
          lut <- setNames(q[[sym_attr]], q$ensembl_gene_id)
          hit_idx <- which(is.na(gene) & !is.na(ids) & ids %in% names(lut))
          
          if (length(hit_idx) > 0) {
            gene[hit_idx] <- unname(lut[ids[hit_idx]])
          }
        }
      }
    }
  }
  
  if (any(is.na(gene)) && .has_gp) {
    q_ids <- unique(ids[is.na(gene) & !is.na(ids)])
    
    if (length(q_ids) > 0) {
      gp <- tryCatch(
        gprofiler2::gconvert(
          q_ids,
          organism = species_organism_code(species),
          target = "SYMBOL",
          mthreshold = Inf,
          filter_na = TRUE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(gp) && nrow(gp) > 0) {
        gp <- gp[!is.na(gp$target) & gp$target != "", , drop = FALSE]
        gp <- gp[!duplicated(gp$input), , drop = FALSE]
        
        lut <- setNames(gp$target, gp$input)
        hit_idx <- which(is.na(gene) & !is.na(ids) & ids %in% names(lut))
        
        if (length(hit_idx) > 0) {
          gene[hit_idx] <- unname(lut[ids[hit_idx]])
        }
      }
    }
  }
  
  gene[is.na(gene) | gene == ""] <- ids[is.na(gene) | gene == ""]
  
  tibble::tibble(
    ID = ids,
    Gene = gene
  )
}

## ============== Backend file helpers (gene lists) ==============
find_backend_file <- function(basename_noext) {
  cands <- c(
    file.path(".",         paste0(basename_noext, ".csv")),
    file.path(".",         paste0(basename_noext, ".tsv")),
    file.path(".",         paste0(basename_noext, ".txt")),
    file.path(".",         paste0(basename_noext, ".xlsx")),
    file.path("/mnt/data", paste0(basename_noext, ".csv")),
    file.path("/mnt/data", paste0(basename_noext, ".tsv")),
    file.path("/mnt/data", paste0(basename_noext, ".txt")),
    file.path("/mnt/data", paste0(basename_noext, ".xlsx"))
  )
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

read_backend_gene_list <- function(basename_noext) {
  f <- find_backend_file(basename_noext)
  if (is.na(f) || is.null(f)) return(character())
  
  ext <- tolower(tools::file_ext(f))
  
  dat <- tryCatch({
    if (ext %in% "csv") {
      readr::read_csv(f, col_names = FALSE, show_col_types = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE)
    } else if (ext %in% "xlsx" && .has_readxl) {
      readxl::read_xlsx(f, col_names = FALSE)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(dat)) return(character())
  unique(na.omit(trimws(as.character(dat[[1]]))))
}

read_backend_gene_table <- function(basename_noext, list_type = c("Symbols", "Ensembl IDs")) {
  list_type <- match.arg(list_type)
  
  f <- find_backend_file(basename_noext)
  if (is.na(f) || !file.exists(f)) {
    return(tibble(GeneRaw = character(), Key = character(), FerroScore = numeric()))
  }
  
  ext <- tolower(tools::file_ext(f))
  
  dat <- tryCatch({
    if (ext %in% "csv") {
      readr::read_csv(f, col_names = TRUE, show_col_types = FALSE)
    } else if (ext %in% c("tsv", "txt")) {
      readr::read_tsv(f, col_names = TRUE, show_col_types = FALSE)
    } else if (ext %in% "xlsx" && .has_readxl) {
      readxl::read_xlsx(f, col_names = TRUE)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(dat)) {
    dat <- tryCatch({
      if (ext %in% "csv") {
        readr::read_csv(f, col_names = FALSE, show_col_types = FALSE)
      } else if (ext %in% c("tsv", "txt")) {
        readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE)
      } else if (ext %in% "xlsx" && .has_readxl) {
        readxl::read_xlsx(f, col_names = FALSE)
      } else {
        NULL
      }
    }, error = function(e) NULL)
  }
  
  if (is.null(dat) || !nrow(dat)) {
    return(tibble(GeneRaw = character(), Key = character(), FerroScore = numeric()))
  }
  
  df <- as.data.frame(dat, check.names = FALSE)
  
  if (ncol(df) == 1) {
    df <- tibble(
      GeneRaw = trimws(as.character(df[[1]])),
      FerroScore = NA_real_
    )
  } else {
    df <- tibble(
      GeneRaw = trimws(as.character(df[[1]])),
      FerroScore = suppressWarnings(as.numeric(df[[2]]))
    )
  }
  
  df$Key <- if (list_type == "Symbols") {
    toupper(df$GeneRaw)
  } else {
    toupper(strip_ensembl_version(df$GeneRaw))
  }
  
  df
}

## ======================= Module CSV helpers ====================
.find_module_file <- function(kind) {
  cands <- c(
    file.path("/mnt/data", paste0(kind, ".csv")),
    file.path(".",         paste0(kind, ".csv"))
  )
  ok <- cands[file.exists(cands)]
  if (length(ok)) ok[1] else NA_character_
}

.read_module_csv <- function(kind) {
  fp <- .find_module_file(kind)
  if (is.na(fp)) return(NULL)
  
  df <- tryCatch(
    readr::read_csv(fp, show_col_types = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(df) || !nrow(df)) return(NULL)
  
  cn <- tolower(names(df))
  if ("module"  %in% cn) names(df)[which(cn == "module")[1]]  <- "Module"
  if ("pathway" %in% cn) names(df)[which(cn == "pathway")[1]] <- "Module"
  if ("gene"    %in% cn) names(df)[which(cn == "gene")[1]]    <- "Gene"
  if (!all(c("Module", "Gene") %in% names(df))) return(NULL)
  
  df %>%
    dplyr::mutate(
      Module = trimws(as.character(Module)),
      Gene   = toupper(trimws(as.character(Gene)))
    ) %>%
    dplyr::filter(nzchar(Module), nzchar(Gene))
}

.modules_from_csv_to_list <- function(kind) {
  df <- .read_module_csv(kind)
  if (is.null(df)) return(NULL)
  split(df$Gene, df$Module) |>
    lapply(function(v) unique(toupper(trimws(v))))
}

## ===================== Enrichment helpers =====================
.read_gmt_T2G <- function(path_no_ext) {
  paths <- c(
    paste0(path_no_ext, ".gmt"),
    file.path("/mnt/data", paste0(basename(path_no_ext), ".gmt"))
  )
  f <- paths[file.exists(paths)][1]
  if (length(f) == 0 || is.na(f)) return(NULL)
  
  gmt <- tryCatch(clusterProfiler::read.gmt(f), error = function(e) NULL)
  if (is.null(gmt) || !nrow(gmt)) return(NULL)
  
  gmt %>%
    dplyr::transmute(
      term = as.character(term),
      gene = toupper(as.character(gene))
    )
}

.csv_modules_to_T2G <- function(dataset_name, modules_list_reactive) {
  lst <- modules_list_reactive()[[dataset_name]]
  if (is.null(lst)) return(NULL)
  
  tibble::tibble(
    term = rep(names(lst), lengths(lst)),
    gene = toupper(unlist(lst, use.names = FALSE))
  )
}

## ==================== NETWORK HELPERS ===================
.get_vst_or_logcpm <- function(counts) {
  if (is.null(counts) || !nrow(counts) || !ncol(counts)) return(NULL)
  
  if (exists("dds", inherits = TRUE)) {
    dds_obj <- tryCatch(get("dds", inherits = TRUE), error = function(e) NULL)
    if (!is.null(dds_obj)) {
      vs <- tryCatch(DESeq2::vst(dds_obj, blind = TRUE), error = function(e) NULL)
      if (!is.null(vs)) {
        return(SummarizedExperiment::assay(vs))
      }
    }
  }
  
  lib <- colSums(counts, na.rm = TRUE)
  lib[is.na(lib) | lib <= 0] <- 1
  cpm <- sweep(counts, 2, lib / 1e6, "/")
  log1p(cpm)
}

.build_corr_edges <- function(vst_mat, genes, cutoff = 0.6, max_edges = 2000) {
  if (is.null(vst_mat) || !nrow(vst_mat) || !ncol(vst_mat)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  gg <- intersect(rownames(vst_mat), toupper(trimws(genes)))
  gg <- gg[!is.na(gg) & gg != ""]
  
  if (length(gg) < 2) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  mat <- vst_mat[gg, , drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  M <- tryCatch(
    stats::cor(t(mat), method = "pearson", use = "pairwise.complete.obs"),
    error = function(e) NULL
  )
  
  if (is.null(M) || !is.matrix(M)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  M[lower.tri(M, diag = TRUE)] <- NA_real_
  
  keep <- which(!is.na(M) & is.finite(M) & abs(M) >= cutoff, arr.ind = TRUE)
  if (!nrow(keep)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  df <- data.frame(
    from   = rownames(M)[keep[, 1]],
    to     = colnames(M)[keep[, 2]],
    weight = as.numeric(M[keep]),
    stringsAsFactors = FALSE
  )
  
  df <- df[is.finite(df$weight), , drop = FALSE]
  df <- df[df$from != df$to, , drop = FALSE]
  
  if (!nrow(df)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  df <- df[order(-abs(df$weight)), , drop = FALSE]
  if (nrow(df) > max_edges) df <- df[seq_len(max_edges), , drop = FALSE]
  
  df
}

.build_string_edges <- function(genes, species = c("Human", "Mouse"),
                                min_score = 400, max_edges = 2000) {
  species <- match.arg(species)
  
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  g <- unique(toupper(trimws(genes)))
  g <- g[!is.na(g) & g != ""]
  if (length(g) < 2) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  sp <- if (species == "Human") 9606L else 10090L
  
  string_db <- tryCatch(
    STRINGdb::STRINGdb$new(
      version = "12",
      species = sp,
      score_threshold = min_score,
      input_directory = ""
    ),
    error = function(e) NULL
  )
  
  if (is.null(string_db)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  df <- data.frame(gene = g, stringsAsFactors = FALSE)
  
  mapped <- tryCatch(
    string_db$map(df, "gene", removeUnmappedRows = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(mapped) || !nrow(mapped) || !"STRING_id" %in% colnames(mapped)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  inter <- tryCatch(
    string_db$get_interactions(mapped$STRING_id),
    error = function(e) NULL
  )
  
  if (is.null(inter) || !nrow(inter)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  lut <- setNames(mapped$gene, mapped$STRING_id)
  inter$from <- lut[inter$from]
  inter$to   <- lut[inter$to]
  
  inter <- inter[
    !is.na(inter$from) & !is.na(inter$to) &
      inter$from != "" & inter$to != "" &
      inter$from != inter$to,
    ,
    drop = FALSE
  ]
  
  if ("combined_score" %in% colnames(inter)) {
    inter <- inter[inter$combined_score >= min_score, , drop = FALSE]
    inter <- unique(inter[, c("from", "to", "combined_score"), drop = FALSE])
    inter <- inter[order(-inter$combined_score), , drop = FALSE]
    if (nrow(inter) > max_edges) inter <- inter[seq_len(max_edges), , drop = FALSE]
    names(inter) <- c("from", "to", "weight")
  } else {
    inter <- unique(inter[, c("from", "to"), drop = FALSE])
    inter$weight <- 1
    if (nrow(inter) > max_edges) inter <- inter[seq_len(max_edges), , drop = FALSE]
  }
  
  inter
}

.build_complete_edges <- function(genes, max_edges = 2000) {
  g <- unique(toupper(genes))
  g <- g[!is.na(g) & g != ""]
  
  if (length(g) < 2) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  
  comb <- utils::combn(g, 2)
  df <- data.frame(
    from = comb[1, ],
    to = comb[2, ],
    weight = 1,
    stringsAsFactors = FALSE
  )
  if (nrow(df) > max_edges) df <- df[seq_len(max_edges), , drop = FALSE]
  df
}

.make_vis_graph <- function(genes, edges_df) {
  nodes <- data.frame(
    id    = unique(toupper(genes)),
    label = unique(toupper(genes)),
    title = unique(toupper(genes)),
    stringsAsFactors = FALSE
  )
  
  edges <- edges_df
  if (is.null(edges) || !nrow(edges)) {
    edges <- data.frame(from = character(), to = character(), weight = numeric())
  }
  
  list(nodes = nodes, edges = edges)
}

## ================== ADD THIS FOR LOGO ACCESS ==================
addResourcePath("assets", normalizePath(".", winslash = "/"))

ui <- fluidPage(
  theme = bslib::bs_theme(
    version     = 5,
    bootswatch  = "flatly",
    primary     = "#0B7285",
    secondary   = "#495057",
    success     = "#2B8A3E",
    info        = "#1971C2",
    warning     = "#E67700",
    danger      = "#C92A2A",
    base_font   = bslib::font_google("Inter"),
    heading_font= bslib::font_google("Plus Jakarta Sans"),
    code_font   = bslib::font_google("JetBrains Mono")
  ),
  
  tags$head(
    tags$style(HTML("
      :root{
        --fe-bg: #f8f9fb;
        --fe-card: #ffffff;
        --fe-border: #e9ecef;
        --fe-muted: #6c757d;
        --fe-ink: #212529;
        --fe-control:  #66c2a5;
        --fe-treated:  #fc8d62;
        --fe-prone:    #d73027;
        --fe-resistant:#2166ac;
      }
      body{ background: var(--fe-bg); color: var(--fe-ink); }

      .page-header h2, .panel-title, .h4, h4 {
        font-weight: 700; letter-spacing: .2px;
      }

      .fe-app-title{
        display:flex;
        align-items:center;
        gap:14px;
        margin-top:6px;
        margin-bottom:6px;
      }

      .fe-app-logo{
        height:72px;
        width:auto;
        display:block;
      }

      .fe-app-title-text{
        font-size:48px;
        font-weight:700;
        line-height:1;
        margin:0;
        color:#212529;
      }

      .fe-card{
        background: var(--fe-card);
        border: 1px solid var(--fe-border);
        border-radius: 14px;
        box-shadow: 0 4px 12px rgba(16,24,40,0.05);
        padding: 18px 18px;
        margin-bottom: 18px;
      }

      .fe-section-title{
        font-weight: 800; font-size: 18px; margin: 0 0 12px 0;
        padding-bottom: 8px; border-bottom: 1px solid var(--fe-border);
      }

      .nav-tabs>li>a{ font-weight:600; color:#495057; }
      .nav-tabs>li.active>a, .nav-tabs>li.active>a:focus, .nav-tabs>li.active>a:hover{
        border-bottom: 3px solid #0B7285; color:#0B7285;
      }

      table.dataTable thead th{ font-weight:700; }
      .dataTables_wrapper .dataTables_filter input{ border-radius:10px; }
      .dataTables_wrapper .dataTables_length select{ border-radius:10px; }

      .btn{ border-radius:10px; font-weight:600; }
      .btn-default{ border-color: var(--fe-border); }

      .tag{
        display:inline-block; padding:2px 8px; border-radius:999px;
        font-size:12px; font-weight:700; color:#fff; margin-right:6px;
      }
      .tag-control{   background: var(--fe-control); }
      .tag-treated{   background: var(--fe-treated); }
      .tag-prone{     background: var(--fe-prone); }
      .tag-resistant{ background: var(--fe-resistant); }

      #legend_panel { 
        border: 1px solid var(--fe-border); 
      }

      /* ===== FIV category legend ===== */
      .fiv-category-box{
        padding-top: 18px;
        padding-left: 10px;
      }

      .fiv-category-box h3{
        font-size: 18px;
        font-weight: 800;
        color: #000000;
        margin-bottom: 18px;
      }

      .fiv-category-row{
        display: flex;
        align-items: center;
        gap: 14px;
        font-size: 18px;
        font-weight: 800;
        color: #000000;
        margin-bottom: 14px;
        line-height: 1.1;
      }

      .fiv-square{
        width: 26px;
        height: 26px;
        display: inline-block;
        border-radius: 6px;
        border: 1px solid rgba(0,0,0,0.15);
        flex-shrink: 0;
      }

      .fiv-square.normal{
        background: #00ff00;
      }

      .fiv-square.mild{
        background: #90ee90;
      }

      .fiv-square.moderate{
        background: #ffff00;
      }

      .fiv-square.high{
        background: #ffa500;
      }

      .fiv-square.severe{
        background: #ff0000;
      }

      /* Old small FIV legend style retained only for backward compatibility */
      .fiv-legend{
        display:flex;
        gap:12px;
        flex-wrap:wrap;
        align-items:center;
        margin-top:10px;
        padding-top:8px;
        border-top:1px solid var(--fe-border);
        font-size:12px;
        color: var(--fe-muted);
        font-weight:700;
      }

      .fiv-item{ 
        display:flex; 
        align-items:center; 
        gap:6px; 
      }

      .fiv-swatch{
        width:14px; 
        height:14px;
        border-radius:4px;
        border:1px solid rgba(0,0,0,0.12);
      }
    "))
  ),
  
  titlePanel(
    div(
      class = "fe-app-title",
      tags$img(
        src = "assets/Logo.png",
        class = "fe-app-logo"
      ),
      tags$div(
        class = "fe-app-title-text",
        "FerroEnrich"
      )
    ),
    windowTitle = "FerroEnrich"
  ),
  
  sidebarLayout(
    sidebarPanel(
      class = "fe-card",
      div(class = "fe-section-title", "Inputs"),
      
      radioButtons(
        "datatype",
        "Data Type",
        c("Read Counts" = "counts"),
        selected = "counts"
      ),
      
      fluidRow(
        column(
          12,
          actionButton(
            "load_demo1",
            "Load Demo Data 1",
            class = "btn-success",
            style = "width:100%; font-weight:bold; margin-bottom:8px;"
          )
        )
      ),
      
      fluidRow(
        column(
          12,
          actionButton(
            "load_demo2",
            "Load Demo Data 2",
            class = "btn-info",
            style = "width:100%; font-weight:bold; margin-bottom:8px;"
          )
        )
      ),
      
      fluidRow(
        column(
          12,
          actionButton(
            "clear_data",
            "Clear Data",
            class = "btn-danger",
            style = "width:100%; font-weight:bold; margin-bottom:12px;"
          )
        )
      ),
      
      textOutput("data_source_text"),
      br(),
      
      fileInput(
        "counts_file",
        "Read counts (genes x samples)",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
      ),
      
      fileInput(
        "meta_file",
        "Sample metadata (rows = samples)",
        accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")
      ),
      
      uiOutput("meta_id_ui"),
      uiOutput("group_col_ui"),
      uiOutput("ref_level_ui"),
      
      selectInput(
        "species",
        "Species (auto-detected, can override)",
        choices = c("Human", "Mouse"),
        selected = "Human"
      ),
      
      tags$small(
        class = "text-muted",
        "Counts: column 1 = Ensembl IDs (ENSG..., ENSMUSG...). ",
        "Metadata: choose which column contains sample IDs."
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        
        tabPanel(
          "Load Data",
          
          div(
            class = "fe-card",
            div(class = "fe-section-title", "About FerroEnrich"),
            
            tags$p(
              strong("FerroEnrich"),
              " is a web-based R Shiny platform designed for integrated analysis of ferroptosis and senescence programs in transcriptomic datasets, with a focus on liver disease biology."
            ),
            
            tags$p(
              "The platform allows users to upload raw RNA-seq count matrices and metadata, perform DESeq2-based differential expression analysis, calculate a Ferroptosis Index Value (FIV), identify ferroptosis-prone and ferroptosis-resistant genes, and explore curated pathway enrichment, interaction networks, and ferroptosis–senescence cross-talk."
            ),
            
            tags$p(
              "FerroEnrich is intended to help researchers connect gene expression changes with biologically meaningful ferroptosis and senescence regulatory modules."
            )
          ),
          
          br(),
          
          div(
            class = "fe-card",
            div(class = "fe-section-title", "Preview"),
            fluidRow(
              column(6, h4("Counts"), DTOutput("counts_head")),
              column(6, h4("Metadata"), DTOutput("meta_head"))
            )
          ),
          
          br(),
          
          div(
            class = "fe-card",
            div(class = "fe-section-title", "Summary"),
            verbatimTextOutput("summary_text")
          ),
          
          br(),
          div(
            class = "fe-card",
            div(class = "fe-section-title", "Contact"),
            
            fluidRow(
              column(
                6,
                h4("Developer"),
                tags$p(strong("Dr. Munichandra Babu Tirumalasetty")),
                tags$p("NYU Grossman School of Medicine"),
                tags$p(
                  "Email: ",
                  tags$a(
                    href = "mailto:munichandra.tirumalasetty@nyulangone.org",
                    "munichandra.tirumalasetty@nyulangone.org"
                  )
                ),
                tags$p(
                  "Alternative Email: ",
                  tags$a(
                    href = "mailto:tmunichandrababu@gmail.com",
                    "tmunichandrababu@gmail.com"
                  )
                )
              ),
              
              column(
                6,
                h4("Principal Investigator"),
                tags$p(strong("Dr. Qing Miao")),
                tags$p("NYU Grossman School of Medicine"),
                tags$p(
                  "Email: ",
                  tags$a(
                    href = "mailto:qing.miao@nyulangone.org",
                    "qing.miao@nyulangone.org"
                  )
                )
              )
            )
          ),
          
          div(
            class = "fe-footer",
            HTML(
              '&copy; 2026 <a href="https://www.miaoqlab.com" target="_blank" rel="noopener noreferrer"><strong>Dr. Miao Lab</strong></a>, NYU Grossman School of Medicine. All rights reserved.'
            )
          )
        ),
          
           # ======================== DEG TAB ========================
        tabPanel(
          "DEG",
          
          div(
            class="fe-card",
            div(class="fe-section-title","DESeq2 Results (design: ~ group)"),
            uiOutput("map_rate"),
            div(
              style = "display:flex; gap:10px; flex-wrap:wrap;",
              actionButton(
                "run_deg",
                "Run DESeq2",
                class = "btn-primary",
                style = "font-weight:bold;"
              ),
              downloadButton("download_deg","Download Results (CSV)")
            ),
            tags$small(
              class = "text-muted",
              "Click Run DESeq2 to start differential expression analysis. A progress popup will appear while running."
            )
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Volcano Settings"),
            fluidRow(
              column(
                3,
                selectInput("sig_metric","Significance metric",
                            choices=c("padj","pvalue"), selected="padj"),
                numericInput("sig_cut","padj/pvalue cutoff",
                             value=0.05, min=1e-10, step=0.01)
              ),
              column(
                3,
                numericInput("lfc_cut","|log2FC| cutoff",
                             value=1, min=0, step=0.1),
                numericInput("top_labels","Label top N genes",
                             value=20, min=0, step=1)
              ),
              column(
                3,
                numericInput("point_size","Point size",
                             value=4.0, min=0.5, step=0.2),
                numericInput("label_size","Label size",
                             value=5.5, min=1, step=0.5)
              ),
              column(
                3,
                numericInput("axis_title_sz","Axis title size", value=20, min=8, step=1),
                numericInput("axis_text_sz","Axis text size",  value=20, min=8, step=1)
              )
            ),
            hr(),
            h4("Custom gene list (optional)"),
            
            fluidRow(
              column(
                6,
                textAreaInput(
                  "gene_text",
                  "Paste gene Symbols",
                  rows = 3,
                  placeholder = "e.g. GPX4, SLC7A11, ACSL4"
                )
              ),
              
              column(
                4,
                br(),
                checkboxInput(
                  "show_matched",
                  "Highlight Prone/Resistant",
                  TRUE
                )
              )
            )
          ),
          
          div(
            class="fe-card",
            plotOutput("volcano", height="580px"),
            div(
              style="margin-top:10px;",
              downloadButton("download_volcano","Download Volcano (PNG)")
            )
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Matched ferroptosis genes"),
            div(
              tags$span(class="tag tag-prone","Prone ↑"),
              tags$span(class="tag tag-resistant","Resistant ↓")
            ),
            uiOutput("ferro_counts")
          ),
          
          div(
            class="fe-card",
            h4("Matched Prone genes"),
            div(style="margin-bottom:8px;",
                downloadButton("dl_prone","Download matched Prone")),
            DTOutput("tbl_prone")
          ),
          
          div(
            class="fe-card",
            h4("Matched Resistant genes"),
            div(style="margin-bottom:8px;",
                downloadButton("dl_resist","Download matched Resistant")),
            DTOutput("tbl_resist")
          ),
          
          div(
            class = "fe-card",
            div(class = "fe-section-title", "Ferroptosis Index Value (FIV)"),
            
            fluidRow(
              column(
                6,
                h4("Classic FIV"),
                plotOutput("fiv_gauge", height = "330px")
              ),
              
              column(
                6,
                br(), br(),
                div(
                  class = "fiv-category-box",
                  h3("FIV category"),
                  
                  div(
                    class = "fiv-category-row",
                    span(class = "fiv-square normal"),
                    span("Normal ≤ 2.50")
                  ),
                  
                  div(
                    class = "fiv-category-row",
                    span(class = "fiv-square mild"),
                    span("Mild 2.50–2.99")
                  ),
                  
                  div(
                    class = "fiv-category-row",
                    span(class = "fiv-square moderate"),
                    span("Moderate 3.00–3.49")
                  ),
                  
                  div(
                    class = "fiv-category-row",
                    span(class = "fiv-square high"),
                    span("High 3.50–4.49")
                  ),
                  
                  div(
                    class = "fiv-category-row",
                    span(class = "fiv-square severe"),
                    span("Severe ≥ 4.50")
                  )
                )
              )
            ),
            
            div(
              style = "display:flex; gap:12px; flex-wrap:wrap; margin-top:10px;",
              downloadButton("dl_fiv_png", "Classic FIV (PNG)"),
              downloadButton("dl_fiv_summary", "Classic FIV (CSV)")
            )
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title",
                "Heatmap of matched Prone ∪ Resistant genes"),
            tags$small(
              class="text-muted",
              HTML(paste0(
                "VST-normalized counts, row z-score, hierarchical clustering. ",
                "<span class='tag tag-prone'>Prone</span>",
                "<span class='tag tag-resistant'>Resistant</span>"
              ))
            ),
            fluidRow(
              column(3, checkboxInput("cluster_rows",  "Cluster rows", TRUE)),
              column(3, checkboxInput("cluster_cols",  "Cluster columns", TRUE)),
              column(3, numericInput("hm_width",  "PNG width (in)", 9, min=6, step=0.5)),
              column(3, numericInput("hm_height", "PNG height (in)", 8, min=6, step=0.5))
            ),
            plotOutput("heatmap", height="700px"),
            div(
              style="margin-top:10px;",
              downloadButton("download_heatmap","Download Heatmap (PNG)")
            )
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","DEG Table"),
            DTOutput("deg_table")
          )
        ),
        
        # ======================== PATHWAYS ========================
        tabPanel(
          "Pathways",
          div(
            class="fe-card",
            div(class="fe-section-title","Module selection"),
            fluidRow(
              column(
                4,
                selectInput(
                  "pathway_sel","Dataset",
                  choices=c("Liver_ferroptosis","Liver_senescence"),
                  selected="Liver_ferroptosis"
                )
              ),
              column(4, uiOutput("module_ui")),
              column(
                2,
                selectInput(
                  "pw_metric","Significance metric",
                  choices=c("padj","pvalue"), selected="padj"
                )
              ),
              column(
                2,
                numericInput("pw_lfc","|log2FC| cutoff", value=0, min=0, step=0.1)
              )
            ),
            fluidRow(
              column(
                3,
                numericInput("pw_cut","padj/pvalue cutoff",
                             value=0.05, min=1e-10, step=0.01)
              ),
              column(
                3,
                actionButton(
                  "run_pathway",
                  "Run Pathway Analysis",
                  class = "btn-primary",
                  style = "width:100%; font-weight:bold; margin-top:25px;"
                )
              )
            )
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Enriched modules"),
            DTOutput("pw_summary"),
            div(
              style="margin-top:6px;",
              downloadButton("dl_pw_summary","Download Summary (CSV)")
            )
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Matched genes in selected module"),
            downloadButton("dl_pw_genes","Download Genes (CSV)"),
            DTOutput("pw_genes")
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Module heatmap"),
            fluidRow(
              column(3, checkboxInput("pw_cluster_rows","Cluster rows", TRUE)),
              column(3, checkboxInput("pw_cluster_cols","Cluster columns", TRUE)),
              column(3, numericInput("pw_hm_w","PNG width (in)", 9, min=6, step=0.5)),
              column(3, numericInput("pw_hm_h","PNG height (in)", 8, min=6, step=0.5))
            ),
            plotOutput("pw_heatmap", height="700px"),
            div(
              style="margin-top:10px;",
              downloadButton("dl_pw_heatmap","Download Module Heatmap (PNG)")
            )
          )
        ),
        
        # ======================= ENRICHMENT =======================
        tabPanel(
          "Enrichment",
          div(
            class="fe-card",
            div(class="fe-section-title","GSEA settings"),
            fluidRow(
              column(
                3,
                selectInput(
                  "enr_dataset","Dataset",
                  choices=c("Liver_ferroptosis","Liver_senescence"),
                  selected="Liver_ferroptosis"
                )
              ),
              column(3, uiOutput("enr_module_ui")),
              column(
                3,
                selectInput(
                  "enr_metric","Ranking metric",
                  choices=c("padj","pvalue"), selected="padj"
                )
              ),
              column(
                3,
                numericInput(
                  "enr_lfc_cut","|log2FC| cutoff (display only)",
                  value=0, min=0, step=0.1
                )
              )
            ),
            fluidRow(
              column(
                3,
                numericInput(
                  "enr_cut","padj/pvalue cutoff (display only)",
                  value=0.05, min=1e-10, step=0.01
                )
              ),
              column(
                3,
                checkboxInput("enr_use_gmt","Use GMT file if available", TRUE)
              ),
              column(
                3,
                actionButton(
                  "run_enr",
                  "Run GSEA",
                  class = "btn-primary",
                  style = "width:100%; font-weight:bold; margin-bottom:8px;"
                ),
                downloadButton("enr_dl_res","Download GSEA results (CSV)")
              )
            ),
            uiOutput("enr_pkg_msg")
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Results"),
            DTOutput("enr_table")
          ),
          
          div(
            class="fe-card",
            div(class="fe-section-title","Enrichment plots"),
            fluidRow(
              column(
                3,
                selectInput("enr_plot_type","Plot type",
                            c("Dot","Bar","Lollipop"),
                            selected="Dot")
              ),
              column(
                3,
                selectInput(
                  "enr_col_by","Color by",
                  c("-log10(padj)","-log10(pvalue)"),
                  selected="-log10(padj)"
                )
              ),
              column(
                3,
                numericInput(
                  "enr_topn","Show top N (by padj)",
                  value=20, min=1, step=1
                )
              ),
              column(
                3,
                checkboxInput("enr_wrap","Wrap module names", TRUE)
              )
            ),
            h5("Upregulated (ES > 0)"),
            plotOutput("enr_plot_up", height="520px"),
            downloadButton("enr_dl_plot_up","Download Upregulated (PNG)"),
            br(), br(),
            h5("Downregulated (ES < 0)"),
            plotOutput("enr_plot_dn", height="520px"),
            downloadButton("enr_dl_plot_dn","Download Downregulated (PNG)")
          )
        ),
        
        # ========================== NETWORK =======================
        tabPanel(
          "Network",
          
          div(
            class = "fe-card",
            div(class = "fe-section-title", "Network"),
            
            fluidRow(
              column(
                3,
                selectInput(
                  "net_dataset",
                  "Gene module dataset",
                  choices = c(
                    "Liver ferroptosis" = "Liver_ferroptosis",
                    "Liver senescence" = "Liver_senescence"
                  ),
                  selected = "Liver_ferroptosis"
                )
              ),
              
              column(
                3,
                uiOutput("net_module_ui")
              ),
              
              column(
                3,
                selectInput(
                  "net_source",
                  "Network source",
                  choices = c(
                    "STRING PPI",
                    "Complete graph"
                  ),
                  selected = "STRING PPI"
                )
              ),
              
              column(
                3,
                numericInput(
                  "net_max_edges",
                  "Max edges",
                  value = 300,
                  min = 50,
                  max = 3000,
                  step = 50
                )
              )
            ),
            
            fluidRow(
              column(
                3,
                numericInput(
                  "net_string_score",
                  "STRING score cutoff",
                  value = 400,
                  min = 150,
                  max = 900,
                  step = 50
                )
              ),
              
              column(
                3,
                br(),
                actionButton("run_network", 
                             "Run Network",
                             class = "btn-primary",
                             style = "width:100%; font-weight:bold; margin-bottom:8px;")
              ),
              
              column(
                3,
                br(),
                actionButton("net_relayout", "Relayout")
              )
            ),
            
            br(),
            
            fluidRow(
              column(
                6,
                downloadButton(
                  "download_network_csv",
                  "Download network table"
                )
              ),
              
              column(
                6,
                downloadButton(
                  "download_current_network_png",
                  "Download current network PNG"
                )
              )
            ),
            
            br(),
            
            div(
              style = "position:absolute; right:35px; z-index:10; background:white; border:1px solid #ddd; border-radius:10px; padding:8px;",
              tags$span(class = "badge bg-danger", "Prone"),
              tags$span(class = "badge bg-primary", "Resistant"),
              uiOutput("net_legend_ui")
            ),
            
            div(
              style = "border:1px solid #e5e7eb; border-radius:18px; padding:10px; background:white;",
              visNetworkOutput("net_vis", height = "900px")
            )
          )
        ),
        # ======================== CROSS-TALK ======================
        tabPanel(
          "Cross-talk",
          
          div(
            class = "fe-card",
            div(class = "fe-section-title", "Cross-talk settings"),
            
            tags$small(
              class = "text-muted",
              HTML(paste0(
                "Visualize cross-talk between <b>ferroptosis modules</b> and ",
                "<b>senescence modules</b> using:<br>",
                "<ul style='padding-left:18px;'>",
                "<li><b>Literature</b> score matrix (−2 .. +2)</li>",
                "<li><b>Correlation</b> of module activity scores</li>",
                "<li><b>Combined</b> = literature × scaled correlation</li>",
                "</ul>"
              ))
            ),
            
            br(), br(),
            
            fluidRow(
              column(
                4,
                radioButtons(
                  "ct_view_type",
                  "Matrix to display",
                  choices = c(
                    "Literature (score_matrix)" = "literature",
                    "Correlation (module-based)" = "correlation",
                    "Combined (Lit × Corr)" = "combined"
                  ),
                  selected = "combined"
                )
              ),
              
              column(
                4,
                selectInput(
                  "ct_order_rows",
                  "Order ferroptosis modules by",
                  choices = c(
                    "Input order",
                    "Hierarchical clustering"
                  ),
                  selected = "Hierarchical clustering"
                )
              ),
              
              column(
                4,
                selectInput(
                  "ct_order_cols",
                  "Order senescence modules by",
                  choices = c(
                    "Input order",
                    "Hierarchical clustering"
                  ),
                  selected = "Hierarchical clustering"
                )
              )
            ),
            
            br(),
            
            fluidRow(
              column(
                12,
                actionButton(
                  "run_crosstalk",
                  "Run Cross-talk Analysis",
                  class = "btn-primary",
                  style = "width:100%; font-weight:bold; margin-top:10px; margin-bottom:8px;"
                ),
                tags$small(
                  class = "text-muted",
                  "Click Run to calculate the cross-talk heatmap, FSI score, top module-pair table, and module network. A progress popup will appear while the analysis is running."
                )
              )
            )
          ),
          
          div(
            class = "fe-card",
            div(
              class = "fe-section-title",
              "Ferroptosis–Senescence module crosstalk heatmap"
            ),
            
            tags$small(
              class = "text-muted",
              HTML(paste0(
                "<b>Red</b> = pro-ferroptotic cross-talk; ",
                "<b>Blue</b> = protective / anti-ferroptotic."
              ))
            ),
            
            br(),
            
            plotOutput("ct_heatmap_matrix", height = "650px"),
            
            div(
              style = "margin-top:10px; display:flex; gap:10px; flex-wrap:wrap;",
              downloadButton("ct_heatmap_matrix_png", "Download heatmap (PNG)"),
              downloadButton("ct_heatmap_matrix_pdf", "Download heatmap (PDF)"),
              downloadButton("ct_matrix_csv", "Download matrix (CSV)")
            )
          ),
          
          div(
            class = "fe-card",
            div(
              class = "fe-section-title",
              "Ferroptosis Susceptibility Index (FSI)"
            ),
            
            tags$small(
              class = "text-muted",
              HTML(paste0(
                "FSI = mean(z-score of pro-ferro modules) − ",
                "mean(z-score of protective modules) per sample."
              ))
            ),
            
            br(),
            
            fluidRow(
              column(
                4,
                uiOutput("fsi_group_col_ui")
              ),
              
              column(
                4,
                checkboxInput(
                  "fsi_show_violin",
                  "Show violin + boxplot when groups are present",
                  TRUE
                )
              ),
              
              column(
                4,
                br(),
                div(
                  style = "display:flex; gap:10px; flex-wrap:wrap;",
                  downloadButton(
                    "fsi_values_csv",
                    "Download FSI values (CSV)"
                  ),
                  downloadButton(
                    "fsi_current_pdf",
                    "Download current FSI PDF"
                  )
                )
              )
            ),
            
            br(),
            
            fluidRow(
              column(
                6,
                h4("FSI distribution"),
                plotOutput("fsi_density", height = "320px")
              ),
              
              column(
                6,
                h4("FSI by group"),
                plotOutput("fsi_by_group", height = "320px")
              )
            )
          ),
          
          div(
            class = "fe-card",
            div(
              class = "fe-section-title",
              "Top ferroptosis–senescence crosstalk pairs"
            ),
            
            tags$small(
              class = "text-muted",
              "Module pairs ranked by |Combined score|."
            ),
            
            br(),
            
            fluidRow(
              column(
                4,
                numericInput(
                  "ct_top_n_edges",
                  "Show top N edges",
                  value = 30,
                  min = 5,
                  step = 5
                )
              ),
              
              column(
                4,
                selectInput(
                  "ct_edge_direction_filter",
                  "Filter by direction",
                  choices = c(
                    "All",
                    "Pro-ferroptotic only",
                    "Protective only"
                  ),
                  selected = "All"
                )
              ),
              
              column(
                4,
                downloadButton(
                  "ct_edges_top_csv",
                  "Download shown edges (CSV)"
                )
              )
            ),
            
            DTOutput("ct_edges_top_tbl")
          ),
          
          div(
            class = "fe-card",
            div(
              class = "fe-section-title",
              "Crosstalk module network"
            ),
            
            tags$small(
              class = "text-muted",
              HTML(paste0(
                "Nodes = modules; edges = top |Combined score| pairs. ",
                "<span class='tag tag-prone'>Red</span> = pro-ferroptotic; ",
                "<span class='tag tag-resistant'>Blue</span> = protective."
              ))
            ),
            
            br(),
            
            fluidRow(
              column(
                4,
                numericInput(
                  "ct_net_top_k",
                  "Use top K edges in network",
                  value = 40,
                  min = 10,
                  step = 5
                )
              ),
              
              column(
                4,
                selectInput(
                  "ct_net_layout",
                  "Network layout",
                  choices = c(
                    "Fruchterman-Reingold" = "fr",
                    "Circle" = "circle",
                    "Kamada-Kawai" = "kk"
                  ),
                  selected = "fr"
                )
              ),
              
              column(
                4,
                downloadButton(
                  "ct_net_pdf_export",
                  "Download network (PDF)"
                )
              )
            ),
            
            plotOutput("ct_net_small", height = "520px")
          )
        ),
        
        # ======================== MANUAL ========================
        tabPanel(
          "Manual",
          div(
            class = "fe-card",
            div(class = "fe-section-title", "FerroEnrich User Manual"),
            tags$p(
              class = "text-muted",
              "The FerroEnrich manual is displayed below. If it does not load, use the button to open it in a new browser tab."
            ),
            tags$a(
              href = "assets/Manual.pdf",
              target = "_blank",
              class = "btn btn-primary",
              "Open Manual PDF"
            ),
            br(), br(),
            tags$iframe(
              src = "assets/Manual.pdf",
              style = "width:100%; height:850px; border:1px solid #e9ecef; border-radius:12px;"
            )
          )
        )
        
      ) # tabsetPanel
    )   # mainPanel
  )     # sidebarLayout
)       # fluidPage

# ---------------- SERVER ----------------
server <- function(input, output, session){
  
  rv <- reactiveValues(
    counts = NULL,
    meta   = NULL,
    source = "No data loaded"
  )
  
  # ---------- Tab click progress popup ----------
  show_tab_progress <- function(tab_name) {
    
    if (!tab_name %in% c("DEG", "Pathways", "Enrichment", "Network", "Cross-talk")) {
      return(NULL)
    }
    
    showNotification(
      paste0("Opening ", tab_name, " tab..."),
      type = "message",
      duration = 2
    )
  }
  
  observeEvent(input$tabs, {
    
    req(input$tabs)
    
    show_tab_progress(input$tabs)
    
  }, ignoreInit = TRUE)
  
  # ---------- File reader ----------
  read_input_file <- function(path, file_name = NULL) {
    
    ext <- tolower(
      tools::file_ext(
        if (is.null(file_name)) path else file_name
      )
    )
    
    tb <- if (ext %in% c("tsv", "txt")) {
      
      readr::read_tsv(
        path,
        col_types = readr::cols()
      )
      
    } else if (ext == "csv") {
      
      readr::read_csv(
        path,
        col_types = readr::cols()
      )
      
    } else if (ext %in% c("xlsx", "xls")) {
      
      validate(
        need(
          .has_readxl,
          "Package 'readxl' is required to read Excel files."
        )
      )
      
      readxl::read_excel(path)
      
    } else {
      
      validate(
        need(
          FALSE,
          "Unsupported file format. Please use CSV, TSV, TXT, XLSX, or XLS."
        )
      )
    }
    
    as.data.frame(tb, check.names = FALSE)
  }
  
  # ---------- helper to find demo files ----------
  find_demo_file <- function(base_name) {
    possible_files <- c(
      file.path(getwd(), paste0(base_name, ".csv")),
      file.path(getwd(), paste0(base_name, ".tsv")),
      file.path(getwd(), paste0(base_name, ".txt")),
      file.path(getwd(), paste0(base_name, ".xlsx")),
      file.path(getwd(), paste0(base_name, ".xls"))
    )
    
    hit <- possible_files[file.exists(possible_files)]
    
    if (length(hit) == 0) {
      stop(paste("Demo file not found in working directory:", base_name))
    }
    
    hit[1]
  }
  
  # ---------- helper to prepare count matrix ----------
  prepare_counts_file <- function(cnt) {
    rn <- cnt[[1]]
    cnt <- cnt[, -1, drop = FALSE]
    names(cnt) <- trimws(names(cnt))
    rownames(cnt) <- trimws(rn)
    cnt
  }
  
  # ---------- Load Demo Data 1 ----------
  observeEvent(input$load_demo1, {
    
    counts_path <- find_demo_file("NASH_1")
    meta_path   <- find_demo_file("Metadata_1")
    
    cnt <- read_input_file(counts_path, basename(counts_path))
    mt  <- read_input_file(meta_path, basename(meta_path))
    
    cnt <- prepare_counts_file(cnt)
    mt[] <- lapply(mt, function(x) if (is.character(x)) trimws(x) else x)
    
    rv$counts <- cnt
    rv$meta   <- mt
    rv$source <- "Demo Data 1 loaded: NASH_1 + Metadata_1"
    
    dds_val(NULL)
    deg_res(NULL)
    map_msg("")
    gc()
    
    updateSelectInput(
      session,
      "species",
      selected = detect_species_from_ids(rownames(cnt)[1])
    )
    
    showNotification("Demo Data 1 loaded successfully.", type = "message")
  })
  
  # ---------- Load Demo Data 2 ----------
  observeEvent(input$load_demo2, {
    
    counts_path <- find_demo_file("NASH_2")
    meta_path   <- find_demo_file("Metadata_2")
    
    cnt <- read_input_file(counts_path, basename(counts_path))
    mt  <- read_input_file(meta_path, basename(meta_path))
    
    cnt <- prepare_counts_file(cnt)
    mt[] <- lapply(mt, function(x) if (is.character(x)) trimws(x) else x)
    
    rv$counts <- cnt
    rv$meta   <- mt
    rv$source <- "Demo Data 2 loaded: NASH_2 + Metadata_2"
    
    dds_val(NULL)
    deg_res(NULL)
    map_msg("")
    gc()
    
    updateSelectInput(
      session,
      "species",
      selected = detect_species_from_ids(rownames(cnt)[1])
    )
    
    showNotification("Demo Data 2 loaded successfully.", type = "message")
  })
  
  # ---------- Clear Data ----------
  observeEvent(input$clear_data, {
    
    rv$counts <- NULL
    rv$meta   <- NULL
    rv$source <- "No data loaded"
    
    dds_val(NULL)
    deg_res(NULL)
    map_msg("")
    
    showNotification(
      "Data cleared. Please upload files or load demo data again.",
      type = "warning"
    )
  })
  
  observeEvent(input$counts_file, {
    req(input$counts_file)
    
    cnt <- read_input_file(input$counts_file$datapath, input$counts_file$name)
    rn <- cnt[[1]]
    cnt <- cnt[, -1, drop = FALSE]
    names(cnt) <- trimws(names(cnt))
    rownames(cnt) <- trimws(rn)
    
    rv$counts <- cnt
    
    if (!is.null(rv$meta)) {
      rv$source <- "User uploaded data"
    } else {
      rv$source <- "Counts uploaded by user"
    }
    
    updateSelectInput(
      session, "species",
      selected = detect_species_from_ids(rownames(cnt)[1])
    )
  })
  
  observeEvent(input$meta_file, {
    req(input$meta_file)
    
    mt <- read_input_file(input$meta_file$datapath, input$meta_file$name)
    mt[] <- lapply(mt, function(x) if (is.character(x)) trimws(x) else x)
    
    rv$meta <- mt
    
    if (!is.null(rv$counts)) {
      rv$source <- "User uploaded data"
    } else {
      rv$source <- "Metadata uploaded by user"
    }
  })
  
  counts_raw <- reactive({
    req(rv$counts)
    rv$counts
  })
  
  meta_raw <- reactive({
    req(rv$meta)
    rv$meta
  })
  
  output$data_source_text <- renderText({
    rv$source
  })
  
  output$meta_id_ui <- renderUI({
    req(meta_raw())
    selectInput("meta_id_col","Sample ID column in metadata",
                choices = names(meta_raw()),
                selected = if ("Description" %in% names(meta_raw())) "Description" else names(meta_raw())[1])
  })
  
  output$group_col_ui <- renderUI({
    req(meta_raw()); choices <- names(meta_raw())
    selectInput("group_col","Group column in metadata",
                choices = choices,
                selected = if ("Type" %in% choices) "Type" else choices[1])
  })
  
  output$ref_level_ui <- renderUI({
    req(meta_raw(), input$group_col)
    lvls <- levels(factor(meta_raw()[[input$group_col]]))
    selectInput("ref_level","Reference level (baseline)", choices = lvls,
                selected = ifelse(any(grepl("^control$", lvls, ignore.case = TRUE)),
                                  lvls[grep("^control$", lvls, ignore.case = TRUE)][1], lvls[1]))
  })
  
  aligned <- reactive({
    req(counts_raw(), meta_raw(), input$meta_id_col, input$group_col)
    cnt  <- counts_raw()
    meta <- meta_raw()
    rownames(meta) <- trimws(meta[[input$meta_id_col]])
    common <- intersect(trimws(colnames(cnt)), rownames(meta))
    validate(need(length(common) >= 2, "No overlapping sample IDs between counts and metadata."))
    cnt  <- cnt[, common, drop = FALSE]
    meta <- meta[common, , drop = FALSE]
    meta$group <- factor(meta[[input$group_col]])
    list(counts = cnt, meta = meta)
  })
  
  output$counts_head <- renderDT({
    df <- counts_raw()
    DT::datatable(head(cbind(ID = rownames(df), df), 10),
                  options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })
  output$meta_head <- renderDT({
    df <- meta_raw()
    DT::datatable(head(df, 10),
                  options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })
  
  output$summary_text <- renderPrint({
    cnt <- aligned()$counts; meta <- aligned()$meta
    cat("Counts:", nrow(cnt), "genes x", ncol(cnt), "samples\n")
    cat("Metadata rows:", nrow(meta), "\n")
    cat("Group column:", input$group_col, "\n")
    cat("Group levels:", paste(levels(meta$group), collapse = ", "), "\n")
    print(table(meta$group))
    cat("\nFirst 5 sample IDs:\n"); print(head(colnames(cnt), 5))
    cat("\nFirst 5 Ensembl IDs:\n"); print(head(rownames(cnt), 5))
  })
  
  # ---------- DEG ----------
  dds_val  <- reactiveVal(NULL)
  deg_res  <- reactiveVal(NULL)
  map_msg  <- reactiveVal("")
  
  observeEvent(input$run_deg, ignoreInit = TRUE, {
    
    withProgress(
      message = "Running DESeq2 differential expression analysis...",
      value = 0,
      {
        
        incProgress(0.05, detail = "Aligning counts and metadata...")
        
        cnt  <- aligned()$counts
        meta <- aligned()$meta
        
        validate(
          need(!is.null(cnt),  "Please load count data first."),
          need(!is.null(meta), "Please load metadata first.")
        )
        
        grp <- meta$group
        
        validate(
          need(length(unique(grp)) >= 2, "At least two groups are required for DESeq2.")
        )
        
        if (!is.null(input$ref_level) && input$ref_level %in% levels(grp)) {
          grp <- stats::relevel(grp, ref = input$ref_level)
        }
        
        meta$group <- grp
        
        incProgress(0.10, detail = "Preparing count matrix...")
        
        cnt_int <- as.matrix(round(cnt))
        mode(cnt_int) <- "integer"
        
        # Replace NA values with zero
        cnt_int[is.na(cnt_int)] <- 0
        
        # Remove genes with all zero counts
        cnt_int <- cnt_int[rowSums(cnt_int) > 0, , drop = FALSE]
        
        incProgress(0.15, detail = "Filtering low-count genes to reduce memory...")
        
        # Memory-friendly low-count filtering for shinyapps.io
        keep <- rowSums(cnt_int >= 10) >= 2
        cnt_int <- cnt_int[keep, , drop = FALSE]
        
        validate(
          need(nrow(cnt_int) >= 100, "Too few genes remain after low-count filtering.")
        )
        
        incProgress(
          0.15,
          detail = paste0("Creating DESeq2 object with ", nrow(cnt_int), " genes...")
        )
        
        dds <- DESeq2::DESeqDataSetFromMatrix(
          countData = cnt_int,
          colData   = meta,
          design    = ~ group
        )
        
        incProgress(
          0.25,
          detail = "Estimating size factors, dispersions, and fitting model..."
        )
        
        dds <- DESeq2::DESeq(
          dds,
          parallel = FALSE,
          quiet = TRUE
        )
        
        incProgress(0.10, detail = "Extracting differential expression results...")
        
        levs <- levels(SummarizedExperiment::colData(dds)$group)
        
        test_level <- if (length(levs) == 2) {
          setdiff(levs, input$ref_level)[1]
        } else {
          levs[length(levs)]
        }
        
        validate(
          need(!is.na(test_level), "Could not identify the comparison group.")
        )
        
        res <- DESeq2::results(
          dds,
          contrast = c("group", test_level, input$ref_level)
        )
        
        res_df <- as.data.frame(res) %>%
          tibble::rownames_to_column("ID") %>%
          dplyr::mutate(ID = strip_ensembl_version(ID))
        
        incProgress(0.15, detail = "Mapping Ensembl IDs to gene symbols...")
        
        species <- isolate(input$species)
        
        map_tbl <- map_ensembl_to_symbol(
          res_df$ID,
          species = species
        )
        
        mapped_n <- sum(!is.na(map_tbl$Gene))
        total_n  <- nrow(map_tbl)
        
        map_msg(
          sprintf(
            "Gene symbol mapping: %d / %d IDs mapped (%.1f%%).",
            mapped_n,
            total_n,
            100 * mapped_n / total_n
          )
        )
        
        incProgress(0.10, detail = "Finalizing DEG table...")
        
        res_out <- res_df %>%
          dplyr::left_join(map_tbl, by = "ID") %>%
          dplyr::relocate(Gene, .after = ID)
        
        dds_val(dds)
        deg_res(res_out)
        
        incProgress(0.05, detail = "DESeq2 analysis completed.")
        
        gc()
      }
    )
  })
  
  output$map_rate <- renderUI({
    req(deg_res())
    
    div(
      style = "margin:6px 0;color:#2c3e50;",
      map_msg()
    )
  })
  
  output$map_rate <- renderUI({
    req(deg_res())
    
    div(
      style = "margin:6px 0;color:#2c3e50;",
      map_msg()
    )
  })
  
  # ---------- custom list (highlight only) ----------
  custom_genes <- reactive({
    
    vals <- character(0)
    
    gene_txt <- input$gene_text
    
    if (!is.null(gene_txt) && length(gene_txt) > 0 && nzchar(trimws(gene_txt))) {
      vals <- unlist(strsplit(gene_txt, "[,;\\s]+"))
    }
    
    vals <- unique(na.omit(trimws(vals)))
    vals <- vals[vals != ""]
    
    vals
  })
  
  # ---------- ferro matching ----------
  # ---------- ferro matching ----------
  ferro_match <- reactive({
    req(deg_res())
    
    res <- deg_res()
    
    sig_res <- res %>%
      dplyr::filter(!is.na(pvalue) & pvalue <= input$sig_cut) %>%
      dplyr::mutate(
        symU = toupper(trimws(Gene)),
        idU  = toupper(strip_ensembl_version(ID))
      )
    
    # Backend list type is fixed to Symbols
    prone_tab  <- read_backend_gene_table("prone_genes", "Symbols")
    resist_tab <- read_backend_gene_table("resistant_genes", "Symbols")
    
    if (nrow(prone_tab) == 0 && nrow(resist_tab) == 0) {
      return(
        list(
          prone = tibble::tibble(),
          resist = tibble::tibble(),
          nprone = 0,
          nresist = 0
        )
      )
    }
    
    m_prone <- sig_res %>%
      dplyr::filter(log2FoldChange > 0) %>%
      dplyr::inner_join(
        prone_tab %>% dplyr::select(Key, FerroScore),
        by = c("symU" = "Key")
      ) %>%
      dplyr::transmute(
        Gene = dplyr::coalesce(Gene, ID),
        log2FoldChange,
        pvalue,
        padj,
        FerroScore
      ) %>%
      dplyr::arrange(dplyr::desc(!is.na(FerroScore)), pvalue)
    
    m_res <- sig_res %>%
      dplyr::filter(log2FoldChange < 0) %>%
      dplyr::inner_join(
        resist_tab %>% dplyr::select(Key, FerroScore),
        by = c("symU" = "Key")
      ) %>%
      dplyr::transmute(
        Gene = dplyr::coalesce(Gene, ID),
        log2FoldChange,
        pvalue,
        padj,
        FerroScore
      ) %>%
      dplyr::arrange(dplyr::desc(!is.na(FerroScore)), pvalue)
    
    list(
      prone = m_prone,
      resist = m_res,
      nprone = nrow(m_prone),
      nresist = nrow(m_res)
    )
  })
  
  output$ferro_counts <- renderUI({
    fm <- ferro_match()
    
    if (fm$nprone == 0 && fm$nresist == 0) {
      div(
        style = "color:#b94a48;",
        "No backend gene lists found or no matches."
      )
    } else {
      div(
        strong("Matched: "),
        span(paste0("Prone = ", fm$nprone, "; Resistant = ", fm$nresist))
      )
    }
  })
  
  output$tbl_prone <- renderDT({
    DT::datatable(
      ferro_match()$prone,
      options = list(scrollX = TRUE, pageLength = 15),
      rownames = FALSE
    )
  })
  
  output$tbl_resist <- renderDT({
    DT::datatable(
      ferro_match()$resist,
      options = list(scrollX = TRUE, pageLength = 15),
      rownames = FALSE
    )
  })
  
  output$dl_prone <- downloadHandler(
    filename = function() {
      "Matched_Prone_genes.csv"
    },
    content = function(file) {
      readr::write_csv(ferro_match()$prone, file)
    }
  )
  
  output$dl_resist <- downloadHandler(
    filename = function() {
      "Matched_Resistant_genes.csv"
    },
    content = function(file) {
      readr::write_csv(ferro_match()$resist, file)
    }
  )
  
  # ----------- Volcano-----------
  # ----------- Volcano -----------
  volcano_data <- reactive({
    req(deg_res())
    
    df <- deg_res()
    
    mcol <- if (identical(input$sig_metric, "padj")) "padj" else "pvalue"
    
    df <- df %>%
      dplyr::mutate(
        neglog10 = -log10(pmax(.data[[mcol]], 1e-300)),
        
        up_sig = 
          (log2FoldChange >= input$lfc_cut) &
          (!is.na(.data[[mcol]]) & .data[[mcol]] <= input$sig_cut),
        
        down_sig = 
          (log2FoldChange <= -input$lfc_cut) &
          (!is.na(.data[[mcol]]) & .data[[mcol]] <= input$sig_cut),
        
        Category = dplyr::case_when(
          up_sig   ~ "Up",
          down_sig ~ "Down",
          TRUE     ~ "NS"
        ),
        
        ID_stripped = strip_ensembl_version(ID),
        symKey      = toupper(trimws(dplyr::coalesce(Gene, ID_stripped))),
        idKey       = toupper(trimws(ID_stripped))
      )
    
    # Custom gene list is now fixed to Symbols only
    glist <- toupper(trimws(custom_genes()))
    glist <- glist[!is.na(glist) & glist != ""]
    
    if (length(glist) > 0) {
      df$Custom <- df$symKey %in% glist
    } else {
      df$Custom <- FALSE
    }
    
    # Matched prone/resistant genes
    fm <- ferro_match()
    
    prone_genes <- toupper(trimws(fm$prone$Gene))
    prone_genes <- prone_genes[!is.na(prone_genes) & prone_genes != ""]
    
    resist_genes <- toupper(trimws(fm$resist$Gene))
    resist_genes <- resist_genes[!is.na(resist_genes) & resist_genes != ""]
    
    df$MatchedProne  <- df$symKey %in% prone_genes
    df$MatchedResist <- df$symKey %in% resist_genes
    
    # Labels
    df$label <- NA_character_
    
    top_n <- suppressWarnings(as.integer(input$top_labels))
    if (is.na(top_n) || top_n < 0) top_n <- 0
    
    if (top_n > 0) {
      topN <- df %>%
        dplyr::filter(Category != "NS") %>%
        dplyr::arrange(dplyr::desc(neglog10)) %>%
        head(top_n)
      
      df$label[match(topN$ID, df$ID)] <- ifelse(
        is.na(topN$Gene) | topN$Gene == "",
        topN$ID,
        topN$Gene
      )
    }
    
    df
  })
  
  
  volc_plot <- reactive({
    df <- volcano_data()
    
    x_finite <- df$log2FoldChange[is.finite(df$log2FoldChange)]
    y_finite <- df$neglog10[is.finite(df$neglog10)]
    
    xmax <- if (length(x_finite)) {
      max(2, ceiling(max(abs(x_finite), na.rm = TRUE)))
    } else {
      2
    }
    
    ymax <- if (length(y_finite)) {
      max(5, ceiling(max(y_finite, na.rm = TRUE)))
    } else {
      5
    }
    
    x_breaks <- pretty(c(-xmax, xmax), n = 7)
    y_breaks <- pretty(c(0, ymax), n = 6)
    
    up_n   <- sum(df$Category == "Up", na.rm = TRUE)
    down_n <- sum(df$Category == "Down", na.rm = TRUE)
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(x = log2FoldChange, y = neglog10)
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = Category),
        alpha = 0.7,
        size = input$point_size
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "Down" = "blue",
          "Up"   = "red",
          "NS"   = "grey75"
        ),
        breaks = c("Up", "Down", "NS")
      ) +
      ggplot2::geom_vline(
        xintercept = c(-input$lfc_cut, input$lfc_cut),
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(input$sig_cut),
        linetype = "dashed",
        color = "black"
      ) +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = label),
        size = input$label_size,
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.25,
        point.padding = 0.15,
        seed = 42,
        na.rm = TRUE
      ) +
      ggplot2::annotate(
        "text",
        x = -xmax * 0.97,
        y = ymax * 0.98,
        label = paste0("\u25BC", down_n),
        color = "blue",
        hjust = 0,
        vjust = 1,
        size = input$label_size + 1,
        fontface = "bold"
      ) +
      ggplot2::annotate(
        "text",
        x = xmax * 0.97,
        y = ymax * 0.98,
        label = paste0("\u25B2", up_n),
        color = "red",
        hjust = 1,
        vjust = 1,
        size = input$label_size + 1,
        fontface = "bold"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(-xmax, xmax),
        breaks = x_breaks,
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, ymax),
        breaks = y_breaks,
        expand = ggplot2::expansion(mult = 0.02)
      ) +
      ggplot2::labs(
        x = "log2 Fold Change",
        y = paste0("-log10(", input$sig_metric, ")"),
        title = "Volcano plot (Up = red, Down = blue)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title = ggplot2::element_text(
          size = input$axis_title_sz,
          face = "bold"
        ),
        axis.text = ggplot2::element_text(
          size = input$axis_text_sz,
          face = "bold"
        ),
        plot.title = ggplot2::element_text(
          size = input$axis_title_sz + 2,
          face = "bold"
        ),
        panel.grid.major = ggplot2::element_line(
          color = "grey85",
          linewidth = 0.4
        ),
        panel.grid.minor = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(
          size = input$axis_text_sz
        )
      )
    
    # Highlight custom pasted Symbols
    if (any(df$Custom, na.rm = TRUE)) {
      p <- p +
        ggplot2::geom_point(
          data = df[df$Custom, ],
          shape = 21,
          stroke = 1.4,
          size = input$point_size + 1.6,
          color = "black",
          fill = NA
        )
    }
    
    # Highlight matched Prone genes
    if (isTRUE(input$show_matched) && any(df$MatchedProne, na.rm = TRUE)) {
      p <- p +
        ggplot2::geom_point(
          data = df[df$MatchedProne, ],
          shape = 24,
          stroke = 1.3,
          size = input$point_size + 1.8,
          color = "#d73027",
          fill = NA
        )
    }
    
    # Highlight matched Resistant genes
    if (isTRUE(input$show_matched) && any(df$MatchedResist, na.rm = TRUE)) {
      p <- p +
        ggplot2::geom_point(
          data = df[df$MatchedResist, ],
          shape = 25,
          stroke = 1.3,
          size = input$point_size + 1.8,
          color = "#2166ac",
          fill = NA
        )
    }
    
    p
  })
  
  
  output$volcano <- renderPlot({
    volc_plot()
  })
  
  
  output$download_volcano <- downloadHandler(
    filename = function() {
      "volcano.png"
    },
    content = function(file) {
      ggplot2::ggsave(
        file,
        plot = volc_plot(),
        width = 8.5,
        height = 6.5,
        dpi = 300
      )
    }
  )
  
  
  filtered_deg <- reactive({
    req(deg_res())
    
    mcol <- if (identical(input$sig_metric, "padj")) "padj" else "pvalue"
    
    deg_res() %>%
      dplyr::mutate(
        sig = 
          (!is.na(.data[[mcol]]) & .data[[mcol]] <= input$sig_cut) &
          (abs(log2FoldChange) >= input$lfc_cut)
      ) %>%
      dplyr::arrange(
        dplyr::desc(sig),
        .data[[mcol]],
        dplyr::desc(abs(log2FoldChange))
      )
  })
  
  
  output$deg_table <- renderDT({
    DT::datatable(
      filtered_deg(),
      options = list(
        scrollX = TRUE,
        pageLength = 25
      ),
      rownames = FALSE
    )
  })
  
  
  output$download_deg <- downloadHandler(
    filename = function() {
      paste0(
        "DESeq2_filtered_",
        input$sig_metric,
        "_",
        input$sig_cut,
        "_lfc",
        input$lfc_cut,
        ".csv"
      )
    },
    content = function(file) {
      readr::write_csv(filtered_deg(), file)
    }
  )
  
  # ---------- Classic FIV (unchanged) ----------
  fiv_metrics <- reactive({
    fm <- ferro_match()
    prone  <- fm$prone  %>% dplyr::mutate(WeightedExpression = log2FoldChange * FerroScore)
    resist <- fm$resist %>% dplyr::mutate(WeightedExpression = log2FoldChange * FerroScore)
    sum_prone     <- sum(prone$WeightedExpression,  na.rm = TRUE)
    sum_resistant <- sum(resist$WeightedExpression, na.rm = TRUE)
    raw_FIV <- sum_prone - sum_resistant
    
    if ((is.na(sum_prone) || sum_prone == 0) && (is.na(sum_resistant) || sum_resistant == 0)) {
      scaled_FIV <- 0
    } else {
      min_FIV <- min(raw_FIV, -sum_resistant, na.rm = TRUE) - 10
      max_FIV <- max(raw_FIV,  sum_prone,     na.rm = TRUE) + 10
      scaled_FIV <- (raw_FIV - min_FIV) / (max_FIV - min_FIV) * 5
      scaled_FIV <- min(max(scaled_FIV, 0), 5)
    }
    
    ferroCategory <- dplyr::case_when(
      scaled_FIV >= 4.50 ~ "Severe",
      scaled_FIV >= 3.50 ~ "High",
      scaled_FIV >= 3.00 ~ "Moderate",
      scaled_FIV >  2.50 ~ "Mild",
      TRUE               ~ "Normal"
    )
    list(sum_prone=sum_prone, sum_resist=sum_resistant, raw_FIV=raw_FIV,
         scaled_FIV=scaled_FIV, category=ferroCategory)
  })
  
  # ---- shared gauge plot helper (used by classic + AI) ----
  .mk_gauge <- function(value_scaled, title_text){
    segs <- tibble::tibble(
      label = c("Normal","Mild","Moderate","High","Severe"),
      xmin  = c(0, 2.5, 3.0, 3.5, 4.5),
      xmax  = c(2.5, 3.0, 3.5, 4.5, 5.0)
    )
    cols <- c(Severe="red", High="orange", Moderate="yellow", Mild="lightgreen", Normal="green")
    ggplot2::ggplot(segs) +
      ggplot2::geom_rect(ggplot2::aes(xmin=xmin, xmax=xmax, ymin=0.55, ymax=1, fill=label),
                         color="white", linewidth=1.1) +
      ggplot2::scale_fill_manual(values=cols, guide="none") +
      ggplot2::scale_x_continuous(limits=c(0,5), expand=c(0,0)) +
      ggplot2::coord_polar(theta="x", start=pi, direction=-1, clip="off") +
      ggplot2::geom_segment(ggplot2::aes(x=value_scaled, xend=value_scaled, y=0.15, yend=1.05),
                            inherit.aes=FALSE, linewidth=1.4, lineend="round", color="black") +
      ggplot2::geom_point(ggplot2::aes(x=value_scaled, y=0.15), inherit.aes=FALSE, size=2, color="black") +
      ggplot2::annotate("text", x=2.5, y=0.02, label=sprintf("FIV\n%.2f", value_scaled),
                        fontface="bold", size=10) +
      ggplot2::theme_void() +
      ggplot2::ggtitle(title_text) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust=0.5, face="bold", size=18),
                     plot.margin = grid::unit(c(2,10,2,10), "pt"))
  }
  
  output$fiv_gauge <- renderPlot({
    fv <- fiv_metrics()
    .mk_gauge(fv$scaled_FIV, paste0("Category: ", fv$category))
  })
  
  output$dl_fiv_summary <- downloadHandler(
    filename = function(){ "FIV_summary.csv" },
    content  = function(file){
      fv <- fiv_metrics()
      out <- tibble::tibble(
        Metric = c("Sum_Prone_Weighted", "Sum_Resistant_Weighted", "Raw_FIV", "Scaled_FIV", "Category"),
        Value  = c(fv$sum_prone, fv$sum_resist, fv$raw_FIV, fv$scaled_FIV, fv$category)
      )
      readr::write_csv(out, file)
    }
  )
  
  # ---- (NEW) high-quality PNG export for classic FIV ----
  output$dl_fiv_png <- downloadHandler(
    filename = function(){ "FIV_gauge.png" },
    content  = function(file){
      fv <- fiv_metrics()
      p <- .mk_gauge(fv$scaled_FIV, paste0("Category: ", fv$category))
      ggplot2::ggsave(file, plot = p, width = 7, height = 5, dpi = 450)
    }
  )
  
  # ---------- (NEW) AI-Enhanced FIV ----------
  # Build light-weight feature frame from DE results + FerroScore matches
  .ai_features <- reactive({
    req(deg_res())
    fm <- ferro_match()
    # core features (robust; no assumptions about python code)
    tibble::tibble(
      n_prone      = nrow(fm$prone),
      n_resist     = nrow(fm$resist),
      sum_prone_ws = sum(fm$prone$log2FoldChange * fm$prone$FerroScore,  na.rm = TRUE),
      sum_res_ws   = sum(fm$resist$log2FoldChange * fm$resist$FerroScore, na.rm = TRUE),
      mean_prone_lfc  = mean(fm$prone$log2FoldChange,  na.rm = TRUE),
      mean_resist_lfc = mean(fm$resist$log2FoldChange, na.rm = TRUE),
      mean_prone_score  = mean(fm$prone$FerroScore,  na.rm = TRUE),
      mean_resist_score = mean(fm$resist$FerroScore, na.rm = TRUE)
    )
  })
  
  # AI prediction (safe, optional)
  ai_fiv <- reactive({
    # disabled or missing python -> return NA and explain
    use_ai <- isTRUE(input$use_ai_fiv %||% TRUE)
    if (!use_ai || is.null(.fa)) {
      return(list(scaled = NA_real_, raw = NA_real_, se = NA_real_, note = "AI disabled or module not found."))
    }
    feats <- .ai_features()
    fv_cl <- fiv_metrics()
    
    # Default prior = classical raw FIV as anchor
    base_raw <- fv_cl$raw_FIV
    
    # Try optional helper(s) if present in ferroai
    # We keep everything wrapped so the app never breaks.
    out <- tryCatch({
      # rank smoothing is optional if user provided function
      if (!is.null(.fa$rank_with_smoothing)) {
        # pass minimal frame: Gene, lfc, pval
        dr <- deg_res() %>% dplyr::transmute(Gene = dplyr::coalesce(Gene, ID),
                                             lfc = log2FoldChange,
                                             p = pvalue)
        invisible(.fa$rank_with_smoothing(reticulate::r_to_py(dr)))
      }
      
      # model contributions (optional)
      cxgb <- if (!is.null(.fa$contrib_xgb))  suppressWarnings( as.numeric(.fa$contrib_xgb(reticulate::r_to_py(feats))) ) else NA_real_
      clog <- if (!is.null(.fa$contrib_logit)) suppressWarnings( as.numeric(.fa$contrib_logit(reticulate::r_to_py(feats))) ) else NA_real_
      
      # Combine: if ensemble available, use it; else blend with classical raw
      if (!is.null(.fa$ensemble_fiv)) {
        res <- .fa$ensemble_fiv(
          base_raw,            # classical raw anchor
          cxgb %||% 0,         # xgb delta
          clog %||% 0,         # logit delta
          as.integer(input$ai_mc %||% 200L)  # MC samples
        )
        # Expect a dict/tuple with fields: raw, scaled, se
        raw    <- as.numeric(res[["raw"]] %||% res[[1]] %||% base_raw)
        scaled <- as.numeric(res[["scaled"]] %||% res[[2]] %||% NA_real_)
        se     <- as.numeric(res[["se"]] %||% res[[3]] %||% NA_real_)
      } else {
        # Simple blend fallback: raw + 0.5*(cxgb + clog)
        raw <- base_raw + 0.5 * ((cxgb %||% 0) + (clog %||% 0))
        # Scale to [0,5] using same scheme as classic
        # guardrails with the same min/max as classic path
        min_FIV <- min(raw, -feats$sum_res_ws, na.rm = TRUE) - 10
        max_FIV <- max(raw,  feats$sum_prone_ws, na.rm = TRUE) + 10
        scaled  <- (raw - min_FIV) / (max_FIV - min_FIV) * 5
        scaled  <- min(max(scaled, 0), 5)
        se <- NA_real_
      }
      
      list(scaled = scaled, raw = raw, se = se, note = "AI model used.")
    }, error = function(e){
      list(scaled = NA_real_, raw = NA_real_, se = NA_real_, note = paste("AI error:", e$message))
    })
    
    out
  })
  
  # Optional text card (can be placed anywhere in UI if you add an output$ai_fiv_card)
  output$ai_fiv_card <- renderUI({
    af <- ai_fiv()
    if (is.na(af$scaled)) return(div(style="color:#777;", "AI-FIV: not available"))
    catg <- dplyr::case_when(
      af$scaled >= 4.50 ~ "Severe",
      af$scaled >= 3.50 ~ "High",
      af$scaled >= 3.00 ~ "Moderate",
      af$scaled >  2.50 ~ "Mild",
      TRUE               ~ "Normal"
    )
    se_txt <- if (is.na(af$se)) "" else sprintf(" (± %.2f SE)", af$se)
    div(
      style="padding:6px 10px; border:1px solid #ddd; border-radius:6px; background:#fafafa;",
      HTML(sprintf("<b>AI-FIV</b>: %.2f%s &nbsp;&nbsp; <b>Category</b>: %s<br/><span style='color:#777'>%s</span>",
                   af$scaled, se_txt, catg, af$note))
    )
  })
  
  # Gauge for AI-FIV (you can add plotOutput('ai_fiv_gauge') in UI if you want a second gauge)
  output$ai_fiv_gauge <- renderPlot({
    af <- ai_fiv()
    req(!is.na(af$scaled))
    .mk_gauge(af$scaled, "AI-enhanced FIV")
  })
  
  # Downloads for AI-FIV
  output$dl_ai_fiv_png <- downloadHandler(
    filename = function(){ "AI_FIV_gauge.png" },
    content  = function(file){
      af <- ai_fiv(); req(!is.na(af$scaled))
      p <- .mk_gauge(af$scaled, "AI-enhanced FIV")
      ggplot2::ggsave(file, plot = p, width = 7, height = 5, dpi = 450)
    }
  )
  output$dl_ai_fiv_csv <- downloadHandler(
    filename = function(){ "AI_FIV_summary.csv" },
    content  = function(file){
      af <- ai_fiv()
      feats <- .ai_features()
      out <- dplyr::bind_cols(
        tibble::tibble(
          Metric = c("AI_raw_FIV", "AI_scaled_FIV", "AI_SE"),
          Value  = c(af$raw, af$scaled, af$se)
        ),
        feats
      )
      readr::write_csv(out, file)
    }
  )
  
  # ---------- Heatmap (matched genes)----------
  make_heatmap <- reactive({
    req(dds_val(), deg_res())
    if (!.has_pheat) return(NULL)
    
    fm <- ferro_match()
    gene_union <- unique(c(fm$prone$Gene, fm$resist$Gene))
    if (length(gene_union) < 2)
      return(list(plot=NULL, message="Not enough matched genes for heatmap."))
    
    dds   <- dds_val()
    meta  <- as.data.frame(SummarizedExperiment::colData(dds))
    
    # Lightweight logCPM instead of re-running VST
    cnt_mat <- DESeq2::counts(dds, normalized = TRUE)
    lib <- colSums(cnt_mat, na.rm = TRUE)
    lib[lib <= 0 | is.na(lib)] <- 1
    cpm <- sweep(cnt_mat, 2, lib / 1e6, "/")
    mat <- log1p(cpm)
    
    map_tbl <- deg_res() %>%
      dplyr::mutate(ID_stripped = strip_ensembl_version(ID),
                    Symbol = dplyr::coalesce(Gene, ID_stripped)) %>%
      dplyr::select(ID_stripped, Symbol)
    
    rn <- strip_ensembl_version(rownames(mat))
    rn_sym <- map_tbl$Symbol[match(rn, map_tbl$ID_stripped)]
    rn_sym[is.na(rn_sym)] <- rn[is.na(rn_sym)]
    rownames(mat) <- rn_sym
    
    keep <- intersect(rownames(mat), gene_union)
    if (length(keep) < 2)
      return(list(plot=NULL, message="Matched genes not present in VST matrix."))
    
    mat <- mat[keep, , drop=FALSE]
    mat_z <- t(scale(t(mat))); mat_z[!is.finite(mat_z)] <- 0
    
    lr <- deg_res() %>%
      dplyr::mutate(Symbol = dplyr::coalesce(Gene, strip_ensembl_version(ID))) %>%
      dplyr::select(Symbol, log2FoldChange)
    lr <- lr[match(rownames(mat_z), lr$Symbol), , drop=FALSE]
    
    set_type <- ifelse(rownames(mat_z) %in% ferro_match()$prone$Gene, "Prone",
                       ifelse(rownames(mat_z) %in% ferro_match()$resist$Gene, "Resistant", ""))
    
    ann_row <- data.frame(Set = set_type,
                          log2FC = lr$log2FoldChange,
                          row.names = rownames(mat_z))
    ann_col <- data.frame(Group = meta$group); rownames(ann_col) <- colnames(mat_z)
    
    lim  <- max(2, stats::quantile(abs(mat_z), 0.98, na.rm=TRUE))
    brks <- seq(-lim, lim, length.out = 101)
    cols <- colorRampPalette(c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C"))(100)
    
    ann_colors <- list(
      Group = c(Control = "#66c2a5", Treated = "#fc8d62"),
      Set   = c(Prone = "#d73027", Resistant = "#2166ac"),
      log2FC = colorRampPalette(c("#2166ac", "#f7f7f7", "#d73027"))(100)
    )
    
    plt <- pheatmap::pheatmap(
      mat_z,
      color  = cols, breaks = brks,
      scale  = "none",
      cluster_rows = isTRUE(input$cluster_rows),
      cluster_cols = isTRUE(input$cluster_cols),
      annotation_row    = ann_row,
      annotation_col    = ann_col,
      annotation_colors = ann_colors,
      annotation_names_row = TRUE,
      annotation_names_col = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize=14, fontsize_row=15, fontsize_col=14, angle_col=90,
      border_color="grey85", na_col="grey95", legend=TRUE,
      treeheight_row=15, treeheight_col=15,
      main = "Matched Prone \u222A Resistant (VST, row z-score)"
    )
    list(plot=plt, message=NULL)
  })
  
  output$heatmap <- renderPlot({
    if (!.has_pheat) { plot.new(); text(0.5,0.5,"Install 'pheatmap' to render heatmaps.", cex=1.3); return() }
    hm <- make_heatmap()
    if (is.null(hm$plot)) { plot.new(); text(0.5,0.5, hm$message %||% "No heatmap to display.", cex=1.2); return() }
    grid::grid.newpage(); grid::grid.draw(hm$plot$gtable)
  })
  
  output$download_heatmap <- downloadHandler(
    filename=function(){ "matched_heatmap.png" },
    content=function(file){
      if (!.has_pheat) return(NULL)
      hm <- make_heatmap()
      if (is.null(hm$plot)) return(NULL)
      png(file, width = input$hm_width, height = input$hm_height, units = "in", res = 300)
      grid::grid.newpage(); grid::grid.draw(hm$plot$gtable)
      dev.off()
    }
  )
  
  
  # ===================== Pathways (same behaviour) =====================
  pathway_modules <- reactive({
    csv_ferr <- .modules_from_csv_to_list("Liver_ferroptosis")
    csv_sene <- .modules_from_csv_to_list("Liver_senescence")
    hard_ferr <- list(
      Hepatic_Ferroptosis_Core = toupper(c(
        "SLC7A11","SLC3A2","GPX4","GCLC","GCLM","GSS","GPX1","PRDX6","TXNRD1","GLS2",
        "SLC1A5","SLC38A1","SLC38A2","SLC25A1","SLC25A11","MDH1","IDH1","ME1")),
      PUFA_Peroxidation_and_Phospholipid_Assembly = toupper(c(
        "ACSL4","ACSL1","ACSL3","ACSL5","LPCAT3","LPCAT1","LPCAT2","MBOAT7","ELOVL2","ELOVL5",
        "FADS1","FADS2","PLA2G4A","ALOX15","ALOX12","ALOX5","PTGS2","POR","CYB5R1","ALOX15B",
        "CYP2C8","CYP2C9","CYP2J2","CYP4F2","CYP4F3","CYP2E1")),
      Ether_Lipid_Peroxisome_Sensitization = toupper(c(
        "AGPS","FAR1","FAR2","GNPAT","TMEM189","PLA2G6","PEX3","PEX5","PEX7","PEX10",
        "PEX11A","PEX11B","PEX16","PEX19","HSD17B4","ACAA1","ACOX1","PECR","ACBD5")),
      Iron_Handling_and_Ferritinophagy = toupper(c(
        "TF","TFRC","TFR2","SLC11A2","SLC40A1","FTH1","FTL","NCOA4","HMOX1","HAMP","CP","STEAP3",
        "SLC39A14","SLC39A8","SFXN1","SFXN2","FXN")),
      CoQ10_FSP1_Axis = toupper(c(
        "AIFM2","NQO1","DHODH","COQ2","COQ3","COQ4","COQ5","COQ6","COQ7","COQ8A","COQ8B",
        "COQ9","COQ10A","COQ10B")),
      GCH1_BH4_Axis = toupper(c("GCH1","QDPR","DHFR","SPR","PTS")),
      NRF2_KEAP1_Antioxidant_Response = toupper(c(
        "NFE2L2","KEAP1","HMOX1","GCLM","GCLC","NQO1","TXNRD1","SOD1","SOD2","PRDX1","PRDX6","HSPB1",
        "TXN","TXN2","TXNIP","GLRX","GLRX2","SRXN1","HMOX2")),
      Mitochondrial_ROS_and_IronSulfur_Protection = toupper(c(
        "CISD1","NFS1","ISCU","FDX1","FDXR","SOD2","PRDX3","GPX4","DLST","NNT","UQCRC2",
        "NDUFS1","NDUFS2","NDUFA9","SDHB","ACO2")),
      p53_SAT1_Lipoxygenase_Link = toupper(c(
        "TP53","CDKN1A","SAT1","ALOX12","ALOX15","GLS2","MDM2","SESN1","SESN2"))
    )
    hard_sene <- list(
      Cell_Cycle_Arrest_p53_p21_p16_pRB = toupper(c(
        "TP53","CDKN1A","CDKN2A","RB1","CDK2","CDK4","CDK6","E2F1","CCND1","CCNE1","CCNA2","MDM2","MDM4")),
      DNA_Damage_Response_DDR_Telomere = toupper(c(
        "ATM","ATR","CHEK1","CHEK2","TP53BP1","H2AFX","PARP1","TERT","TERF1","TERF2","POT1","RTEL1","DCLRE1C","WRN")),
      SASP_Cytokines_Chemokines = toupper(c(
        "IL6","CXCL8","IL1A","IL1B","TNF","CCL2","CCL5","CXCL1","CXCL2","CXCL3","CXCL10","LIF","CSF2")),
      SASP_GrowthFactors_TGF_Signals = toupper(c(
        "TGFB1","TGFB2","INHBA","VEGFA","FGF2","HGF","EGF","AREG","EREG","IGFBP3","IGFBP7")),
      SASP_Proteases_ECM_Remodeling = toupper(c(
        "MMP1","MMP3","MMP9","MMP12","MMP14","TIMP1","TIMP2","SERPINE1","PLAU","COL1A1","COL3A1","FN1")),
      Mitochondrial_Dysfunction_ROS = toupper(c(
        "SOD2","PRDX3","PRDX5","GPX1","UQCRC2","NDUFS1","NDUFA9","MT-ND1","MT-CO1","TFAM","PPARGC1A","SIRT3")),
      Autophagy_Lysosome = toupper(c(
        "BECN1","ATG5","ATG7","ATG3","ATG12","ATG16L1","ULK1","MAP1LC3B","SQSTM1","LAMP1","LAMP2","CTSB","CTSL","CTSD","LYZ")),
      NFkB_STING_Pathway_Regulators = toupper(c(
        "RELA","NFKB1","NFKBIA","IKBKB","TAB1","TAB2","TAB3","MAP3K7","MB21D1","TMEM173","TBK1","IRF3")),
      Immune_Clearance_of_Senescent_Cells = toupper(c(
        "MICA","MICB","ULBP1","ULBP2","ULBP3","KLRK1","CX3CL1","ICAM1","VCAM1","GZMB","PRF1","IFNG")),
      Hepatic_Stellate_Cell_Senescence_and_Fibrosis = toupper(c(
        "PDGFRB","PDGFRL","ACTA2","TAGLN","COL1A1","COL1A2","COL3A1","LAMA2","LAMB1","LOX","LOXL2","YAP1","TEAD1","CTGF","TGFBR1","TGFBR2"))
    )
    list(
      Liver_ferroptosis = csv_ferr %||% hard_ferr,
      Liver_senescence  = csv_sene %||% hard_sene
    )
  })
  
  pathway_sets <- reactive({
    pm <- pathway_modules()
    lapply(pm, function(lst) unique(unlist(lst, use.names = FALSE)))
  })
  
  output$module_ui <- renderUI({
    mods <- names(pathway_modules()[[ input$pathway_sel ]] %||% list())
    selectInput("module_sel","Module", choices = c("All", mods), selected = "All")
  })
  observeEvent(input$pathway_sel, {
    mods <- names(pathway_modules()[[ input$pathway_sel ]] %||% list())
    updateSelectInput(session, "module_sel", choices = c("All", mods), selected = "All")
  })
  
  summarize_pathway <- function(res_tbl, genes_uc, metric="padj", cut=0.05, lfc_cut=1){
    mcol <- if (metric=="padj") "padj" else "pvalue"
    tbl  <- res_tbl %>% dplyr::mutate(Symbol = toupper(dplyr::coalesce(Gene, strip_ensembl_version(ID))))
    in_set <- tbl$Symbol %in% genes_uc
    is_sig <- (!is.na(tbl[[mcol]]) & tbl[[mcol]] <= cut & abs(tbl$log2FoldChange) >= lfc_cut)
    up     <- is_sig & tbl$log2FoldChange > 0
    down   <- is_sig & tbl$log2FoldChange < 0
    mk_fisher <- function(flag){
      a <- sum(in_set & flag); b <- sum(!in_set & flag)
      c <- sum(in_set & !flag); d <- sum(!in_set & !flag)
      stats::fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater")$p.value
    }
    p_up   <- mk_fisher(up)
    p_down <- mk_fisher(down)
    list(
      summary = tibble::tibble(
        Direction = c("Down regulated","Up regulated"),
        Pval      = c(p_down, p_up),
        nGenes    = c(sum(in_set & down), sum(in_set & up))
      ),
      genes = tbl %>% dplyr::filter(in_set) %>%
        dplyr::transmute(Gene = Symbol, log2FoldChange, pvalue, padj)
    )
  }
  
  pw_current <- eventReactive(input$run_pathway, {
    
    withProgress(
      message = "Running pathway module analysis...",
      value = 0,
      {
        
        incProgress(0.05, detail = "Checking DEG results...")
        
        validate(
          need(!is.null(deg_res()), "Please run DESeq2 first from the DEG tab.")
        )
     
        res_tbl <- deg_res()
        
        incProgress(0.10, detail = "Reading selected pathway settings...")
        
        metric <- isolate(input$pw_metric)
        cut    <- isolate(input$pw_cut)
        lfc    <- isolate(input$pw_lfc)
        ds     <- isolate(input$pathway_sel)
        
        incProgress(0.10, detail = "Loading pathway modules...")
        
        mods    <- pathway_modules()[[ds]]
        all_set <- pathway_sets()[[ds]]
        
        validate(
          need(!is.null(mods) && length(mods) > 0,
               "No pathway modules found for the selected dataset.")
        )
        
        rows <- list()
        genes_by_mod <- list()
        
        module_names <- names(mods)
        n_mods <- length(module_names)
        
        incProgress(0.10, detail = paste0("Testing ", n_mods, " modules..."))
        
        for (i in seq_along(module_names)) {
          
          nm <- module_names[i]
          
          out <- summarize_pathway(
            res_tbl,
            mods[[nm]],
            metric,
            cut,
            lfc
          )
          
          rows[[nm]] <- out$summary %>%
            dplyr::mutate(Module = nm)
          
          genes_by_mod[[nm]] <- out$genes %>%
            dplyr::mutate(Module = nm)
          
          incProgress(
            0.50 / max(n_mods, 1),
            detail = paste0("Completed module ", i, " of ", n_mods, ": ", nm)
          )
        }
        
        incProgress(0.10, detail = "Summarizing all modules...")
        
        out_all <- summarize_pathway(
          res_tbl,
          all_set,
          metric,
          cut,
          lfc
        )
        
        rows[["All"]] <- out_all$summary %>%
          dplyr::mutate(Module = "All")
        
        genes_by_mod[["All"]] <- out_all$genes %>%
          dplyr::mutate(Module = "All")
        
        incProgress(0.10, detail = "Adjusting p-values...")
        
        all_sum <- dplyr::bind_rows(rows)
        all_sum$adj.Pval <- p.adjust(all_sum$Pval, method = "BH")
        
        incProgress(0.05, detail = "Finalizing pathway results...")
        
        pathway_output <- list(
          sum = all_sum %>%
            dplyr::select(Direction, adj.Pval, nGenes, Module),
          genes = dplyr::bind_rows(genes_by_mod)
        )
        
        gc()
        
        incProgress(0.05, detail = "Pathway analysis completed.")
        
        pathway_output
      }
    )
    
  }, ignoreInit = TRUE)
  
  output$pw_summary <- renderDT({
    req(input$run_pathway > 0)
    req(pw_current())
    df <- pw_current()$sum %>% dplyr::filter(Module == input$module_sel)
    DT::datatable(df, options=list(scrollX=TRUE, pageLength=20), rownames=FALSE)
  })
  output$dl_pw_summary <- downloadHandler(
    filename=function(){ paste0("module_summary_", input$pathway_sel, "_", input$module_sel, ".csv") },
    content=function(file){ readr::write_csv(pw_current()$sum %>% dplyr::filter(Module==input$module_sel), file) }
  )
  output$pw_genes <- renderDT({
    req(pw_current())
    g <- pw_current()$genes %>% dplyr::filter(Module == input$module_sel)
    DT::datatable(g, options=list(scrollX=TRUE, pageLength=25), rownames=FALSE)
  })
  output$dl_pw_genes <- downloadHandler(
    filename=function(){ paste0("module_genes_", input$pathway_sel, "_", input$module_sel, ".csv") },
    content=function(file){ readr::write_csv(pw_current()$genes %>% dplyr::filter(Module==input$module_sel), file) }
  )
  
  make_pw_heatmap <- reactive({
    req(dds_val(), deg_res()); if (!.has_pheat) return(NULL)
    ds  <- input$pathway_sel
    mod <- input$module_sel
    genes <- if (mod=="All") pathway_sets()[[ds]] else pathway_modules()[[ds]][[mod]]
    if (length(genes) < 2) return(list(plot=NULL, message="Too few genes in set."))
    dds   <- dds_val()
    meta  <- as.data.frame(SummarizedExperiment::colData(dds))
    vst_m <- DESeq2::vst(dds, blind = TRUE)
    mat   <- SummarizedExperiment::assay(vst_m)
    map_tbl <- deg_res() %>%
      dplyr::mutate(ID_stripped = strip_ensembl_version(ID),
                    Symbol = toupper(dplyr::coalesce(Gene, ID_stripped))) %>%
      dplyr::select(ID_stripped, Symbol)
    rn <- strip_ensembl_version(rownames(mat))
    rn_sym <- toupper(map_tbl$Symbol[match(rn, map_tbl$ID_stripped)])
    rn_sym[is.na(rn_sym)] <- toupper(rn[is.na(rn_sym)])
    rownames(mat) <- rn_sym
    keep <- intersect(rownames(mat), genes)
    if (length(keep) < 2) return(list(plot=NULL, message="Selected module genes not found in VST matrix."))
    mat <- mat[keep, , drop=FALSE]
    mat_z <- t(scale(t(mat))); mat_z[!is.finite(mat_z)] <- 0
    lr <- deg_res() %>%
      dplyr::mutate(Symbol = toupper(dplyr::coalesce(Gene, strip_ensembl_version(ID)))) %>%
      dplyr::select(Symbol, log2FoldChange)
    lr <- lr[match(rownames(mat_z), lr$Symbol), , drop=FALSE]
    ann_row <- data.frame(log2FC = lr$log2FoldChange, row.names = rownames(mat_z))
    ann_col <- data.frame(Group = meta$group); rownames(ann_col) <- colnames(mat_z)
    lim  <- max(2, stats::quantile(abs(mat_z), 0.98, na.rm=TRUE))
    brks <- seq(-lim, lim, length.out = 101)
    cols <- colorRampPalette(c("#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61","#D7191C"))(100)
    ann_colors <- list(
      Group = c(Control = "#66c2a5", Treated = "#fc8d62"),
      log2FC = colorRampPalette(c("#2166ac", "#f7f7f7", "#d73027"))(100)
    )
    plt <- pheatmap::pheatmap(
      mat_z, color=cols, breaks=brks, scale="none",
      cluster_rows = isTRUE(input$pw_cluster_rows),
      cluster_cols = isTRUE(input$pw_cluster_cols),
      annotation_row=ann_row, annotation_col=ann_col,
      annotation_colors=ann_colors,
      show_rownames=TRUE, show_colnames=TRUE,
      fontsize=14, fontsize_row=15, fontsize_col=14, angle_col=90,
      border_color="grey85", na_col="grey95",
      main=paste0("Module: ", input$module_sel, " (VST, row z-score)")
    )
    list(plot=plt, message=NULL)
  })
  
  output$pw_heatmap <- renderPlot({
    if (!.has_pheat) { plot.new(); text(0.5,0.5,"Install 'pheatmap' to render heatmaps.", cex=1.3); return() }
    hm <- make_pw_heatmap()
    if (is.null(hm$plot)) { plot.new(); text(0.5,0.5, hm$message %||% "No heatmap to display.", cex=1.2); return() }
    grid::grid.newpage(); grid::grid.draw(hm$plot$gtable)
  })
  
  output$dl_pw_heatmap <- downloadHandler(
    filename=function(){ paste0("heatmap_", input$pathway_sel, "_", input$module_sel, ".png") },
    content=function(file){
      if (!.has_pheat) return(NULL)
      hm <- make_pw_heatmap(); if (is.null(hm$plot)) return(NULL)
      png(file, width = input$pw_hm_w, height = input$pw_hm_h, units = "in", res = 300)
      grid::grid.newpage(); grid::grid.draw(hm$plot$gtable)
      dev.off()
    }
  )
  
  # ====================== ENRICHMENT (GSEA) ======================
  
  # ---------- helpers for GMT/CSV term2gene ----------
  .read_gmt_T2G <- function(path_no_ext){
    paths <- c(paste0(path_no_ext, ".gmt"),
               file.path("/mnt/data", paste0(basename(path_no_ext), ".gmt")))
    f <- paths[file.exists(paths)][1]
    if (length(f) == 0 || is.na(f)) return(NULL)
    gmt <- tryCatch(clusterProfiler::read.gmt(f), error=function(e) NULL)
    if (is.null(gmt) || !nrow(gmt)) return(NULL)
    gmt %>% dplyr::transmute(term = as.character(term),
                             gene = toupper(as.character(gene)))
  }
  
  .csv_modules_to_T2G <- function(dataset_name, modules_list_reactive){
    lst <- modules_list_reactive()[[dataset_name]]
    if (is.null(lst)) return(NULL)
    tibble::tibble(
      term = rep(names(lst), lengths(lst)),
      gene = toupper(unlist(lst, use.names = FALSE))
    )
  }
  
  # ---------- GSEA inputs ----------
  enr_T2G <- reactive({
    ds <- input$enr_dataset %||% "Liver_ferroptosis"
    if (isTRUE(input$enr_use_gmt)) {
      t2g <- .read_gmt_T2G(file.path("/mnt/data", ds))
      if (!is.null(t2g)) return(t2g)
    }
    .csv_modules_to_T2G(ds, pathway_modules)
  })
  
  output$enr_module_ui <- renderUI({
    t2g <- isolate(enr_T2G())
    mods <- if (is.null(t2g)) character(0) else unique(t2g$term)
    selectInput("enr_module","Module", choices = c("All", mods), selected = "All")
  })
  
  output$enr_pkg_msg <- renderUI({
    if (!requireNamespace("clusterProfiler", quietly = TRUE))
      return(div(style="color:#b94a48;",
                 "Packages 'clusterProfiler' and 'enrichplot' are required for GSEA."))
    HTML("")
  })
  
  # ---------- Build ranked gene list from DE results ----------
  .build_geneList <- function(res_tbl, metric = c("padj","pvalue")){
    metric <- match.arg(metric)
    mcol   <- if (metric == "padj") "padj" else "pvalue"
    df <- res_tbl %>%
      dplyr::mutate(
        ID_stripped = strip_ensembl_version(ID),
        SymbolRaw   = dplyr::coalesce(trimws(Gene), ID_stripped),
        Symbol      = toupper(trimws(SymbolRaw))
      ) %>%
      dplyr::filter(!is.na(.data[[mcol]]),
                    is.finite(log2FoldChange),
                    !is.na(Symbol), nzchar(Symbol)) %>%
      dplyr::distinct(Symbol, .keep_all = TRUE)
    
    stat <- if ("stat" %in% names(df) && all(is.finite(df$stat))) df$stat
    else sign(df$log2FoldChange) * -log10(pmax(df[[mcol]], 1e-300))
    
    rank_tbl <- tibble::tibble(Symbol = df$Symbol, stat = stat) %>%
      dplyr::filter(!is.na(stat), is.finite(stat), !is.na(Symbol), nzchar(Symbol)) %>%
      dplyr::group_by(Symbol) %>%
      dplyr::summarise(stat = stat[which.max(abs(stat))], .groups = "drop")
    
    gl <- stats::setNames(rank_tbl$stat, rank_tbl$Symbol)
    gl <- gl[!is.na(names(gl)) & nzchar(names(gl))]
    validate(need(length(gl) >= 10, "Too few ranked genes after cleanup for GSEA."))
    sort(gl, decreasing = TRUE)
  }
  
  # ---------- Run GSEA ----------
  enr_run <- eventReactive(input$run_enr, {
    
    withProgress(
      message = "Running GSEA enrichment analysis...",
      value = 0,
      {
        
        incProgress(0.05, detail = "Checking DEG results...")
        
        validate(
          need(!is.null(deg_res()), "Please run DESeq2 first from the DEG tab.")
        )
      
        validate(
          need(requireNamespace("clusterProfiler", quietly = TRUE),
               "Package 'clusterProfiler' is required for GSEA.")
        )
        
        incProgress(0.10, detail = "Preparing TERM2GENE module list...")
        
        t2g <- isolate(enr_T2G())
        
        validate(
          need(!is.null(t2g) && nrow(t2g) > 0,
               "No module definitions found. Please check GMT/CSV module files.")
        )
        
        t2g <- t2g %>%
          dplyr::mutate(
            gene = toupper(trimws(gene)),
            term = trimws(term)
          ) %>%
          dplyr::filter(
            nzchar(gene),
            nzchar(term),
            !is.na(gene),
            !is.na(term)
          )
        
        incProgress(0.15, detail = "Building ranked gene list from DESeq2 results...")
        
        gl <- .build_geneList(
          deg_res(),
          metric = isolate(input$enr_metric)
        )
        
        incProgress(0.15, detail = "Matching ranked genes with module genes...")
        
        t2g <- t2g %>%
          dplyr::filter(gene %in% names(gl))
        
        validate(
          need(nrow(t2g) > 0,
               "No module genes overlap with the ranked gene list.")
        )
        
        incProgress(0.15, detail = "Running GSEA calculation...")
        
        set.seed(123)
        
        gsea_res <- clusterProfiler::GSEA(
          geneList     = gl,
          TERM2GENE    = t2g,
          minGSSize    = 5,
          maxGSSize    = 300,
          pvalueCutoff = 1.0,
          verbose      = FALSE,
          by           = "fgsea"
        )
        
        incProgress(0.25, detail = "Finalizing GSEA results...")
        
        gc()
        
        incProgress(0.15, detail = "GSEA completed.")
        
        gsea_res
      }
    )
    
  }, ignoreInit = TRUE)
  # ---------- Convert GSEA result to table ----------
  enr_tbl <- reactive({
    g  <- enr_run()
    df <- as.data.frame(g)
    if (!nrow(df)) return(df)
    
    have <- names(df)
    qval <- if ("qvalues" %in% have) df[["qvalues"]]
    else if (requireNamespace("qvalue", quietly = TRUE)) {
      tryCatch(qvalue::qvalue(df$pvalue)$qvalues, error = function(e) rep(NA_real_, nrow(df)))
    } else rep(NA_real_, nrow(df))
    
    res <- tibble::tibble(
      Module          = df$ID,
      Description     = df$Description,
      setSize         = df$setSize,
      enrichmentScore = df$enrichmentScore,
      NES             = df$NES,
      pvalue          = df$pvalue,
      padj            = df$p.adjust,
      qvalue          = qval,
      rank            = if ("rank" %in% have) df$rank else NA_integer_,
      leading_edge    = if ("leading_edge" %in% have) df$leading_edge else NA_character_,
      core_enrichment = if ("core_enrichment" %in% have) df$core_enrichment else NA_character_,
      Dataset         = input$enr_dataset
    )
    
    selected_module <- isolate(input$enr_module)
    
    if (!is.null(selected_module) && selected_module != "All")
      res <- dplyr::filter(res, Module == selected_module)
    
    dplyr::arrange(res, padj, dplyr::desc(NES))
  })
  
  output$enr_table <- renderDT({
    req(input$run_enr > 0)
    req(enr_tbl())
    
    DT::datatable(
      enr_tbl(),
      options = list(
        scrollX = TRUE,
        pageLength = 15,
        order = list(list(6, 'asc'), list(4, 'desc'))
      ),
      rownames = FALSE
    )
  })
  output$enr_dl_res <- downloadHandler(
    filename = function(){ paste0("GSEA_", input$enr_dataset, "_results_full.csv") },
    content  = function(file){ readr::write_csv(enr_tbl(), file) }
  )
  
  # ===================== PLOTS: Up & Down panels =====================
  
  # Safe wrapping (no simplify=)
  .enr_wrap_labels <- function(x, do_wrap, width = 28){
    if (isTRUE(do_wrap)) stringr::str_wrap(x, width = width) else x
  }
  
  # Build plot-ready df, then split by NES sign
  .enr_plot_df_both <- reactive({
    df <- enr_tbl(); req(!is.null(df), nrow(df) > 0)
    
    # core/leading-edge counts
    nCore <- ifelse(!is.na(df$core_enrichment) & nzchar(df$core_enrichment),
                    vapply(strsplit(df$core_enrichment, "/"), length, integer(1)),
                    0L)
    
    # color metric selection
    col_choice <- input$enr_col_by %||% "-log10(padj)"
    col_by <- if (grepl("pvalue", col_choice, ignore.case = TRUE)) "pvalue" else "padj"
    pvec   <- if (identical(col_by, "padj")) df$padj else df$pvalue
    col_log <- -log10(pmax(pvec, 1e-300))
    color_lab <- if (identical(col_by,"padj")) "-log10(padj)" else "-log10(pvalue)"
    
    # base tidy frame
    tidy <- tibble::tibble(
      Module      = df$Module,
      Description = df$Description,
      nCore       = nCore,
      setSize     = df$setSize,
      NES         = df$NES,
      padj        = df$padj,
      pvalue      = df$pvalue,
      color_val   = col_log
    )
    
    # split
    up <- tidy %>% dplyr::filter(NES >  0)
    dn <- tidy %>% dplyr::filter(NES <  0)
    
    # order & topN by padj then NES within each subset
    topN <- max(1L, suppressWarnings(as.integer(input$enr_topn %||% 20L)))
    keep_top <- function(d){
      if (!nrow(d)) return(d)
      ord <- order(d$padj, -d$NES, na.last = TRUE)
      d[ head(ord, min(length(ord), topN)), , drop=FALSE]
    }
    up <- keep_top(up); dn <- keep_top(dn)
    
    # factor order (top->bottom) and wrap labels
    relevel_wrap <- function(d){
      if (!nrow(d)) return(d)
      d$Description <- factor(d$Description, levels = rev(d$Description))
      lab <- .enr_wrap_labels(levels(d$Description), input$enr_wrap)
      d$Description <- factor(as.character(d$Description),
                              levels = levels(d$Description), labels = lab)
      d
    }
    up <- relevel_wrap(up); dn <- relevel_wrap(dn)
    
    list(up = up, dn = dn, color_lab = color_lab)
  })
  
  # Common plot builder (type: Dot / Bar / Lollipop)
  .enr_build_plot <- function(df, type, color_lab, title){
    # theme
    base_theme <- ggplot2::theme_minimal(base_size = 18) +
      ggplot2::theme(
        axis.text       = ggplot2::element_text(color = "black", face = "bold"),
        axis.title      = ggplot2::element_text(color = "black", face = "bold"),
        axis.text.x     = ggplot2::element_text(size = 16, color = "black", face = "bold"),
        axis.text.y     = ggplot2::element_text(size = 16, color = "black", face = "bold"),
        plot.title      = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
        legend.title    = ggplot2::element_text(face = "bold", color = "black"),
        legend.text     = ggplot2::element_text(color = "black"),
        panel.border    = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.2),
        legend.position = "right"
      )
    
    col_grad  <- ggplot2::scale_color_gradient(name = color_lab, low = "blue", high = "red")
    fill_grad <- ggplot2::scale_fill_gradient (name = color_lab, low = "blue", high = "red")
    
    if (!nrow(df)) {
      return(ggplot2::ggplot() + ggplot2::theme_void() +
               ggplot2::annotate("text", x=0.5, y=0.5, label="No terms to display",
                                 size=6, fontface="bold") +
               base_theme)
    }
    
    if (identical(type, "Dot")) {
      ggplot2::ggplot(df, ggplot2::aes(x = nCore, y = Description)) +
        ggplot2::geom_point(ggplot2::aes(size = nCore, color = color_val)) +
        ggplot2::scale_size_continuous(name = "N of Genes", range = c(3, 10)) +
        col_grad +
        ggplot2::labs(x = "N of Genes", y = NULL, title = title) +
        base_theme
      
    } else if (identical(type, "Lollipop")) {
      ggplot2::ggplot(df, ggplot2::aes(y = Description, x = nCore)) +
        ggplot2::geom_segment(ggplot2::aes(yend = Description, x = 0, xend = nCore),
                              linewidth = 1.1, color = "grey60") +
        ggplot2::geom_point(ggplot2::aes(color = color_val), size = 4) +
        col_grad +
        ggplot2::labs(x = "N of Genes", y = NULL, title = title) +
        base_theme
      
    } else { # Bar
      ggplot2::ggplot(df, ggplot2::aes(x = nCore, y = Description, fill = color_val)) +
        ggplot2::geom_col(width = 0.7) +
        fill_grad +
        ggplot2::labs(x = "N of Genes", y = NULL, title = title) +
        base_theme
    }
  }
  
  # ---- Render up/down plots ----
  output$enr_plot_up <- renderPlot({
    req(input$run_enr > 0)
    d <- .enr_plot_df_both()
    .enr_build_plot(
      d$up,
      input$enr_plot_type %||% "Dot",
      d$color_lab,
      "Upregulated (ES > 0)"
    )
  })
  
  output$enr_plot_dn <- renderPlot({
    req(input$run_enr > 0)
    d <- .enr_plot_df_both()
    .enr_build_plot(
      d$dn,
      input$enr_plot_type %||% "Dot",
      d$color_lab,
      "Downregulated (ES < 0)"
    )
  })
  
  # ---- Downloads ----
  output$enr_dl_plot_up <- downloadHandler(
    filename = function(){ paste0("enrichment_up_", input$enr_plot_type, ".png") },
    content  = function(file){
      d <- .enr_plot_df_both()
      p <- .enr_build_plot(d$up, input$enr_plot_type %||% "Dot",
                           d$color_lab, "Upregulated (ES > 0)")
      ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
  
  output$enr_dl_plot_dn <- downloadHandler(
    filename = function(){ paste0("enrichment_down_", input$enr_plot_type, ".png") },
    content  = function(file){
      d <- .enr_plot_df_both()
      p <- .enr_build_plot(d$dn, input$enr_plot_type %||% "Dot",
                           d$color_lab, "Downregulated (ES < 0)")
      ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
  ## ===================== NETWORK: helpers + server =====================
  
  # Palette for log2FC
  .lfc_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))
  
  
  # Pull gene and log2FC from DEG table
  deg_map <- reactive({
    
    df <- tryCatch(deg_res(), error = function(e) NULL)
    
    if (is.null(df) || !nrow(df)) {
      return(tibble::tibble(gene = character(), lfc = numeric()))
    }
    
    nm <- tolower(names(df))
    
    gene_vec <- if ("gene" %in% nm) {
      df[[which(nm == "gene")[1]]]
    } else if ("symbol" %in% nm) {
      df[[which(nm == "symbol")[1]]]
    } else {
      rownames(df)
    }
    
    lfc_col <- if ("log2foldchange" %in% nm) {
      names(df)[which(nm == "log2foldchange")[1]]
    } else if ("log2fc" %in% nm) {
      names(df)[which(nm == "log2fc")[1]]
    } else if ("lfc" %in% nm) {
      names(df)[which(nm == "lfc")[1]]
    } else {
      NULL
    }
    
    if (is.null(lfc_col)) {
      return(tibble::tibble(gene = character(), lfc = numeric()))
    }
    
    tibble::tibble(
      gene = toupper(as.character(gene_vec)),
      lfc  = as.numeric(df[[lfc_col]])
    ) |>
      dplyr::distinct(gene, .keep_all = TRUE)
  })
  
  
  # Complete graph fallback
  .build_complete_edges <- function(genes, max_edges = 2000) {
    
    g <- unique(toupper(genes))
    
    if (length(g) < 2) {
      return(data.frame(
        from = character(),
        to = character(),
        weight = numeric()
      ))
    }
    
    comb <- utils::combn(g, 2)
    
    df <- data.frame(
      from = comb[1, ],
      to = comb[2, ],
      weight = 1,
      stringsAsFactors = FALSE
    )
    
    if (nrow(df) > max_edges) {
      df <- df[seq_len(max_edges), , drop = FALSE]
    }
    
    df
  }
  
  
  # STRING edges
  .build_string_edges <- function(genes,
                                  species = c("Human", "Mouse"),
                                  min_score = 400,
                                  max_edges = 2000) {
    
    species <- match.arg(species)
    
    if (!requireNamespace("STRINGdb", quietly = TRUE)) {
      return(NULL)
    }
    
    sp <- if (species == "Human") 9606L else 10090L
    
    sdb <- tryCatch(
      STRINGdb::STRINGdb$new(
        version = "12",
        species = sp,
        score_threshold = min_score,
        input_directory = ""
      ),
      error = function(e) NULL
    )
    
    if (is.null(sdb)) return(NULL)
    
    df <- data.frame(
      gene = unique(toupper(genes)),
      stringsAsFactors = FALSE
    )
    
    mapped <- tryCatch(
      sdb$map(df, "gene", removeUnmappedRows = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(mapped) || !nrow(mapped)) return(NULL)
    
    ints <- tryCatch(
      sdb$get_interactions(mapped$STRING_id),
      error = function(e) NULL
    )
    
    if (is.null(ints) || !nrow(ints)) return(NULL)
    
    lut <- setNames(mapped$gene, mapped$STRING_id)
    
    ints$from <- lut[ints$from]
    ints$to   <- lut[ints$to]
    
    ints <- ints[!is.na(ints$from) & !is.na(ints$to), , drop = FALSE]
    ints <- unique(ints[c("from", "to", "combined_score")])
    
    if (!nrow(ints)) return(NULL)
    
    ints <- ints[order(-ints$combined_score), , drop = FALSE]
    
    if (nrow(ints) > max_edges) {
      ints <- ints[seq_len(max_edges), , drop = FALSE]
    }
    
    names(ints) <- c("from", "to", "weight")
    
    ints
  }
  
  
  # Nodes and edges shell
  .make_vis_graph <- function(genes, edges_df) {
    
    genes <- unique(toupper(genes))
    
    nodes <- data.frame(
      id = genes,
      label = genes,
      title = genes,
      stringsAsFactors = FALSE
    )
    
    edges <- edges_df
    
    if (is.null(edges) || !nrow(edges)) {
      edges <- data.frame(
        from = character(),
        to = character(),
        weight = numeric()
      )
    }
    
    list(nodes = nodes, edges = edges)
  }
  
  
  ## ----------------- NETWORK MODULE SERVER -----------------
  
  .has_vis    <- requireNamespace("visNetwork", quietly = TRUE)
  .has_igraph <- requireNamespace("igraph", quietly = TRUE)
  
  
  output$net_module_ui <- renderUI({
    
    ds <- input$net_dataset %||% "Liver_ferroptosis"
    mods <- names(pathway_modules()[[ds]] %||% list())
    
    selectInput(
      "net_module",
      "Module",
      choices = c("All", mods),
      selected = "All"
    )
  })
  
  
  # Seed for reproducible layout
  .net_seed <- reactiveVal(123L)
  
  observeEvent(input$net_relayout, {
    .net_seed(sample.int(1e6, 1))
  })
  
  
  # Reactive gene list
  net_genes <- reactive({
    
    ds <- input$net_dataset %||% "Liver_ferroptosis"
    mods <- pathway_modules()[[ds]]
    
    validate(
      need(!is.null(mods), "No modules available for the selected dataset.")
    )
    
    if (is.null(input$net_module) || input$net_module == "All") {
      unique(toupper(unlist(mods, use.names = FALSE)))
    } else {
      unique(toupper(mods[[input$net_module]]))
    }
  })
  
  
  # Build edges: STRING PPI or Complete graph only
  net_edges <- eventReactive(input$run_network, {
    
    withProgress(
      message = "Building network...",
      value = 0,
      {
        incProgress(0.05, detail = "Preparing module gene list...")
        
        validate(
          need(!is.null(deg_res()), "Please run DESeq2 first from the DEG tab.")
        )
        
        gs <- isolate(net_genes())
        
        validate(
          need(length(gs) >= 2, "Too few genes to build a network. Need at least 2 genes.")
        )
        
        src  <- isolate(input$net_source %||% "STRING PPI")
        maxE <- as.integer(isolate(input$net_max_edges %||% 300))
        
        incProgress(0.10, detail = paste0("Selected source: ", src))
        
        # ---------------- STRING PPI ----------------
        if (identical(src, "STRING PPI")) {
          
          sc <- as.integer(isolate(input$net_string_score %||% 400))
          
          incProgress(
            0.30,
            detail = paste0("Retrieving STRING interactions with score ≥ ", sc, "...")
          )
          
          ed <- tryCatch(
            .build_string_edges(
              genes = gs,
              species = isolate(input$species %||% "Human"),
              min_score = sc,
              max_edges = maxE
            ),
            error = function(e) NULL
          )
          
          if (!is.null(ed) && nrow(ed) > 0) {
            incProgress(0.50, detail = "STRING network completed.")
            gc()
            return(ed)
          }
          
          incProgress(0.20, detail = "No STRING edges found. Using fallback complete graph.")
          
          ed <- .build_complete_edges(gs, max_edges = maxE)
          
          incProgress(0.50, detail = "Fallback network completed.")
          gc()
          
          return(ed)
        }
        
        # ---------------- COMPLETE GRAPH ----------------
        incProgress(0.40, detail = "Building complete graph network...")
        
        ed <- .build_complete_edges(gs, max_edges = maxE)
        
        incProgress(0.50, detail = "Complete graph network completed.")
        gc()
        
        ed
      }
    )
    
  }, ignoreInit = TRUE)
  
  
  # Combine nodes and edges
  net_graph <- reactive({
    
    req(input$run_network > 0)
    
    gs <- isolate(net_genes())
    ed <- net_edges()
    
    dat <- .make_vis_graph(gs, ed)
    
    nodes <- dat$nodes
    edges <- dat$edges
    
    nodes$id <- make.unique(as.character(nodes$id))
    
    edges <- edges[
      edges$from %in% nodes$id & edges$to %in% nodes$id,
      ,
      drop = FALSE
    ]
    
    list(nodes = nodes, edges = edges)
  })
  
  
  # Color nodes by log2FC
  node_colors_by_lfc <- reactive({
    
    dat <- net_graph()
    nodes <- dat$nodes
    df <- deg_map()
    
    m <- dplyr::left_join(
      tibble::tibble(
        id = nodes$id,
        label_uc = toupper(nodes$label)
      ),
      df |>
        dplyr::rename(label_uc = gene),
      by = "label_uc"
    )
    
    lfc <- suppressWarnings(as.numeric(m$lfc))
    n <- nrow(nodes)
    
    if (n == 0L || all(is.na(lfc))) {
      return(list(
        bg = rep("#E3F2FD", n),
        vmax = 0,
        pal = NULL,
        has = FALSE
      ))
    }
    
    vmax <- max(abs(lfc), na.rm = TRUE)
    vmax <- ifelse(vmax > 0, vmax, 1)
    
    s <- (lfc / vmax + 1) / 2
    
    pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
    idx <- pmax(1, pmin(101, round(s * 100) + 1))
    
    cols <- pal[idx]
    cols[is.na(cols)] <- "#E3F2FD"
    
    list(
      bg = cols,
      vmax = vmax,
      pal = pal,
      has = TRUE
    )
  })
  
  
  # Compact legend
  output$net_legend_ui <- renderUI({
    
    colinfo <- node_colors_by_lfc()
    
    if (!isTRUE(colinfo$has)) return(NULL)
    
    tags$div(
      tags$strong("log₂FC scale"),
      br(),
      tags$span(style = sprintf("color:%s;", colinfo$pal[1]), "↓ Low  "),
      "|",
      tags$span(style = sprintf("color:%s;", colinfo$pal[101]), "  High ↑"),
      br(),
      sprintf("±%.2f", colinfo$vmax)
    )
  })
  
  
  # Interactive visNetwork only
  output$net_vis <- visNetwork::renderVisNetwork({
    
    req(input$run_network > 0)
    
    dat <- net_graph()
    nodes <- dat$nodes
    edges <- dat$edges
    colinfo <- node_colors_by_lfc()
    
    vn_nodes <- data.frame(
      id = nodes$id,
      label = nodes$label,
      title = paste0("<b>", nodes$label, "</b>"),
      color = colinfo$bg,
      value = 1,
      borderWidth = 1,
      shadow = TRUE,
      stringsAsFactors = FALSE
    )
    
    vn_edges <- if (nrow(edges)) {
      data.frame(
        from = edges$from,
        to = edges$to,
        value = abs(edges$weight %||% 1),
        color = "#9E9E9E",
        smooth = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        from = character(),
        to = character()
      )
    }
    
    visNetwork::visNetwork(
      vn_nodes,
      vn_edges,
      height = "900px",
      width = "100%"
    ) %>%
      visNetwork::visNodes(
        shape = "dot",
        font = list(
          size = 28,
          face = "bold",
          color = "black"
        ),
        scaling = list(min = 12, max = 28)
      ) %>%
      visNetwork::visEdges(
        arrows = list(to = FALSE),
        smooth = FALSE
      ) %>%
      visNetwork::visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = TRUE
      ) %>%
      visNetwork::visIgraphLayout(
        randomSeed = .net_seed()
      ) %>%
      visNetwork::visPhysics(
        stabilization = TRUE
      ) %>%
      visNetwork::visLayout(
        randomSeed = .net_seed()
      ) %>%
      visNetwork::visInteraction(
        navigationButtons = TRUE,
        zoomView = TRUE,
        dragView = TRUE
      ) %>%
      visNetwork::visExport(
        type = "png",
        name = "FerroEnrich_network",
        label = "Save network PNG"
      )
  })
  
  
  # Download network edge table
  # Download network edge table
  output$download_network_csv <- downloadHandler(
    
    filename = function() {
      src <- input$net_source %||% "Network"
      src <- gsub("[^A-Za-z0-9]+", "_", src)
      paste0("FerroEnrich_", src, "_network_edges_", Sys.Date(), ".csv")
    },
    
    content = function(file) {
      
      req(input$run_network > 0)
      
      dat <- net_graph()
      edges <- dat$edges
      
      if (is.null(edges) || !nrow(edges)) {
        edges <- data.frame(
          from = character(),
          to = character(),
          weight = numeric()
        )
      }
      
      write.csv(edges, file, row.names = FALSE)
    }
  )
  
  
  # Download CURRENT selected network as high-quality PNG
  output$download_current_network_png <- downloadHandler(
    
    filename = function() {
      src <- input$net_source %||% "Network"
      src <- gsub("[^A-Za-z0-9]+", "_", src)
      paste0("FerroEnrich_", src, "_network_", Sys.Date(), ".png")
    },
    
    content = function(file) {
      
      req(input$run_network > 0)
      req(.has_igraph)
      
      dat <- net_graph()
      nodes <- dat$nodes
      edges <- dat$edges
      
      validate(
        need(!is.null(nodes) && nrow(nodes) > 0, "No network nodes available.")
      )
      
      nodes$id <- as.character(nodes$id)
      nodes$label <- as.character(nodes$label)
      
      if (is.null(edges) || !nrow(edges)) {
        edges <- data.frame(
          from = character(),
          to = character(),
          weight = numeric(),
          stringsAsFactors = FALSE
        )
      }
      
      edges$from <- as.character(edges$from)
      edges$to   <- as.character(edges$to)
      
      # Remove invalid edges and self-loops
      edges <- edges[
        edges$from %in% nodes$id &
          edges$to %in% nodes$id &
          edges$from != edges$to,
        ,
        drop = FALSE
      ]
      
      # Build igraph object
      if (nrow(edges) > 0) {
        
        g <- igraph::graph_from_data_frame(
          d = edges[, c("from", "to"), drop = FALSE],
          directed = FALSE,
          vertices = nodes
        )
        
      } else {
        
        g <- igraph::make_empty_graph(n = nrow(nodes), directed = FALSE)
        igraph::V(g)$name <- nodes$id
      }
      
      # Stable layout matching current network seed
      set.seed(.net_seed())
      
      if (igraph::ecount(g) > 0) {
        lay <- igraph::layout_with_fr(g)
      } else {
        lay <- igraph::layout_in_circle(g)
      }
      
      # Node colors by log2FC
      df_lfc <- deg_map()
      node_names <- igraph::V(g)$name
      
      lfc_vec <- df_lfc$lfc[
        match(
          toupper(node_names),
          toupper(df_lfc$gene)
        )
      ]
      
      if (all(is.na(lfc_vec))) {
        
        node_cols <- rep("#E3F2FD", length(node_names))
        vmax <- NA_real_
        
      } else {
        
        vmax <- max(abs(lfc_vec), na.rm = TRUE)
        vmax <- ifelse(is.finite(vmax) && vmax > 0, vmax, 1)
        
        pal <- colorRampPalette(
          c("#2166AC", "#F7F7F7", "#B2182B")
        )(101)
        
        scaled <- (lfc_vec / vmax + 1) / 2
        idx <- pmax(1, pmin(101, round(scaled * 100) + 1))
        
        node_cols <- pal[idx]
        node_cols[is.na(node_cols)] <- "#E3F2FD"
      }
      
      # Edge width scaling
      if (nrow(edges) > 0 && "weight" %in% names(edges)) {
        
        edge_w <- suppressWarnings(as.numeric(edges$weight))
        edge_w[is.na(edge_w)] <- 1
        
        max_w <- max(abs(edge_w), na.rm = TRUE)
        max_w <- ifelse(is.finite(max_w) && max_w > 0, max_w, 1)
        
        edge_w_scaled <- 1 + 4 * abs(edge_w) / max_w
        edge_w_scaled[!is.finite(edge_w_scaled)] <- 1
        
      } else {
        
        edge_w_scaled <- 1
      }
      
      # Size settings for large vs small networks
      n_nodes <- igraph::vcount(g)
      
      node_size <- if (n_nodes > 150) {
        8
      } else if (n_nodes > 100) {
        10
      } else if (n_nodes > 70) {
        13
      } else {
        18
      }
      
      label_size <- if (n_nodes > 150) {
        0.65
      } else if (n_nodes > 100) {
        0.80
      } else if (n_nodes > 70) {
        1.00
      } else {
        1.25
      }
      
      title_txt <- paste0(
        "FerroEnrich ",
        input$net_source %||% "Network",
        " network"
      )
      
      # High-quality PNG output
      grDevices::png(
        filename = file,
        width = 16,
        height = 14,
        units = "in",
        res = 600
      )
      
      par(
        mar = c(1, 1, 4, 1),
        xpd = TRUE
      )
      
      plot(
        g,
        layout = lay,
        vertex.color = node_cols,
        vertex.frame.color = "grey35",
        vertex.size = node_size,
        vertex.label = igraph::V(g)$name,
        vertex.label.cex = label_size,
        vertex.label.color = "black",
        edge.color = "grey55",
        edge.width = edge_w_scaled,
        main = title_txt
      )
      
      if (is.finite(vmax)) {
        legend(
          "bottomright",
          legend = c(
            "Low log2FC",
            "High log2FC",
            paste0("Scale ±", round(vmax, 2))
          ),
          pch = c(21, 21, NA),
          pt.bg = c("#2166AC", "#B2182B", NA),
          pt.cex = c(2, 2, NA),
          bty = "n",
          cex = 0.9
        )
      }
      
      grDevices::dev.off()
    }
  )
  
  ## ===================== CROSS-TALK (module ↔ module) =====================
 
    ct_ferro_mods <- c(
    "Hepatic_Ferroptosis_Core",
    "PUFA_Peroxidation_and_Phospholipid_Assembly",
    "Ether_Lipid_Peroxisome_Sensitization",
    "Iron_Handling_and_Ferritinophagy",
    "CoQ10_FSP1_Axis",
    "GCH1_BH4_Axis",
    "NRF2_KEAP1_Antioxidant_Response",
    "Mitochondrial_ROS_and_IronSulfur_Protection",
    "p53_SAT1_Lipoxygenase_Link"
  )
  
  ct_sene_mods <- c(
    "Cell_Cycle_Arrest_p53_p21_p16_pRB",
    "DNA_Damage_Response_DDR_Telomere",
    "SASP_Cytokines_Chemokines",
    "SASP_GrowthFactors_TGF_Signals",
    "SASP_Proteases_ECM_Remodeling",
    "Mitochondrial_Dysfunction_ROS",
    "Autophagy_Lysosome",
    "NFkB_STING_Pathway_Regulators",
    "Immune_Clearance_of_Senescent_Cells",
    "Hepatic_Stellate_Cell_Senescence_and_Fibrosis"
  )
  
  
  # ---- literature-based score matrix (-2 .. +2) ----
  ct_literature_mat <- local({
    
    m <- matrix(
      0,
      nrow = length(ct_ferro_mods),
      ncol = length(ct_sene_mods),
      dimnames = list(ct_ferro_mods, ct_sene_mods)
    )
    
    ferro_protective <- c(
      "Hepatic_Ferroptosis_Core",
      "CoQ10_FSP1_Axis",
      "GCH1_BH4_Axis",
      "NRF2_KEAP1_Antioxidant_Response",
      "Mitochondrial_ROS_and_IronSulfur_Protection"
    )
    
    ferro_proferro <- setdiff(ct_ferro_mods, ferro_protective)
    
    sene_proferro <- c(
      "Cell_Cycle_Arrest_p53_p21_p16_pRB",
      "DNA_Damage_Response_DDR_Telomere",
      "SASP_Cytokines_Chemokines",
      "SASP_GrowthFactors_TGF_Signals",
      "SASP_Proteases_ECM_Remodeling",
      "Mitochondrial_Dysfunction_ROS",
      "NFkB_STING_Pathway_Regulators",
      "Hepatic_Stellate_Cell_Senescence_and_Fibrosis"
    )
    
    sene_protective <- c(
      "Autophagy_Lysosome",
      "Immune_Clearance_of_Senescent_Cells"
    )
    
    for (f in ct_ferro_mods) {
      for (s in ct_sene_mods) {
        if (f %in% ferro_proferro && s %in% sene_proferro) {
          m[f, s] <- 2
        } else if (f %in% ferro_proferro && s %in% sene_protective) {
          m[f, s] <- 1
        } else if (f %in% ferro_protective && s %in% sene_proferro) {
          m[f, s] <- -2
        } else if (f %in% ferro_protective && s %in% sene_protective) {
          m[f, s] <- -1
        } else {
          m[f, s] <- 0
        }
      }
    }
    
    m["p53_SAT1_Lipoxygenase_Link", "Cell_Cycle_Arrest_p53_p21_p16_pRB"] <- 2
    m["p53_SAT1_Lipoxygenase_Link", "DNA_Damage_Response_DDR_Telomere"] <- 2
    m["PUFA_Peroxidation_and_Phospholipid_Assembly", "SASP_Cytokines_Chemokines"] <- 2
    m["PUFA_Peroxidation_and_Phospholipid_Assembly", "Mitochondrial_Dysfunction_ROS"] <- 2
    m["Iron_Handling_and_Ferritinophagy", "Autophagy_Lysosome"] <- 1
    
    m
  })
  
  
  # ---- module activity scores per sample ----
  ct_module_scores <- eventReactive(input$run_crosstalk, {
    
    withProgress(
      message = "Running ferroptosis-senescence cross-talk analysis...",
      value = 0,
      {
        incProgress(0.05, detail = "Checking input data...")
        req(aligned())
        
        incProgress(0.10, detail = "Loading ferroptosis and senescence modules...")
        
        pm <- pathway_modules()
        ferro_list <- pm$Liver_ferroptosis
        sene_list  <- pm$Liver_senescence
        
        validate(
          need(
            !is.null(ferro_list) && !is.null(sene_list),
            "No ferroptosis/senescence module definitions available."
          )
        )
        
        incProgress(0.15, detail = "Preparing normalized expression matrix...")
        
        counts <- aligned()$counts
        logcpm <- .get_vst_or_logcpm(as.matrix(counts))
        
        validate(
          need(!is.null(logcpm), "Could not compute expression matrix for cross-talk analysis.")
        )
        
        sp <- isolate(input$species %||% "Human")
        
        incProgress(0.15, detail = "Mapping gene IDs to symbols...")
        
        map_tbl <- map_ensembl_to_symbol(
          rownames(logcpm),
          species = sp
        )
        
        gene_sym <- toupper(map_tbl$Gene)
        nSamp <- ncol(logcpm)
        
        validate(
          need(nSamp >= 2, "At least two samples are required for cross-talk correlation.")
        )
        
        score_mat <- function(mod_names, mod_list, module_type = "module") {
          
          out <- matrix(
            NA_real_,
            nrow = length(mod_names),
            ncol = nSamp,
            dimnames = list(mod_names, colnames(logcpm))
          )
          
          for (i in seq_along(mod_names)) {
            
            mname <- mod_names[i]
            gset  <- toupper(mod_list[[mname]] %||% character(0))
            
            idx <- which(
              !is.na(gene_sym) &
                gene_sym %in% gset
            )
            
            if (length(idx) >= 2) {
              out[i, ] <- colMeans(
                logcpm[idx, , drop = FALSE],
                na.rm = TRUE
              )
            }
            
            incProgress(
              0.20 / max(length(mod_names), 1),
              detail = paste0("Scoring ", module_type, " module: ", mname)
            )
          }
          
          out
        }
        
        incProgress(0.10, detail = "Scoring ferroptosis modules...")
        
        ferro_scores <- score_mat(
          ct_ferro_mods,
          ferro_list,
          module_type = "ferroptosis"
        )
        
        incProgress(0.10, detail = "Scoring senescence modules...")
        
        sene_scores <- score_mat(
          ct_sene_mods,
          sene_list,
          module_type = "senescence"
        )
        
        incProgress(0.15, detail = "Finalizing module scores...")
        
        out <- list(
          ferro = ferro_scores,
          sene = sene_scores,
          samples = colnames(logcpm)
        )
        
        gc()
        
        incProgress(0.10, detail = "Cross-talk module scoring completed.")
        
        out
      }
    )
    
  }, ignoreInit = TRUE)
  
  
  # ---- correlation matrix between ferro and senescence modules ----
  ct_corr_mat <- reactive({
    
    req(input$run_crosstalk > 0)
    
    withProgress(
      message = "Computing cross-talk correlation matrix...",
      value = 0,
      {
        incProgress(0.15, detail = "Reading module scores...")
        
        sc <- ct_module_scores()
        
        Fm <- sc$ferro
        Sm <- sc$sene
        
        if (
          is.null(Fm) || is.null(Sm) ||
          ncol(Fm) < 2 || ncol(Sm) < 2
        ) {
          return(
            matrix(
              NA_real_,
              nrow = length(ct_ferro_mods),
              ncol = length(ct_sene_mods),
              dimnames = list(ct_ferro_mods, ct_sene_mods)
            )
          )
        }
        
        mat <- matrix(
          NA_real_,
          nrow = nrow(Fm),
          ncol = nrow(Sm),
          dimnames = list(rownames(Fm), rownames(Sm))
        )
        
        total_steps <- nrow(Fm) * nrow(Sm)
        step_count <- 0
        
        incProgress(0.10, detail = "Calculating module-pair correlations...")
        
        for (i in seq_len(nrow(Fm))) {
          for (j in seq_len(nrow(Sm))) {
            
            mat[i, j] <- suppressWarnings(
              stats::cor(
                Fm[i, ],
                Sm[j, ],
                use = "pairwise.complete.obs"
              )
            )
            
            step_count <- step_count + 1
            
            incProgress(
              0.60 / max(total_steps, 1),
              detail = paste0("Correlation pair ", step_count, " of ", total_steps)
            )
          }
        }
        
        incProgress(0.15, detail = "Correlation matrix completed.")
        
        gc()
        
        mat
      }
    )
  })
  
  
  # ---- bundle all matrices ----
  ct_mats_all <- reactive({
    
    req(input$run_crosstalk > 0)
    
    lit <- ct_literature_mat
    cor <- ct_corr_mat()
    
    if (!all(dim(lit) == dim(cor))) {
      cor <- matrix(
        NA_real_,
        nrow = nrow(lit),
        ncol = ncol(lit),
        dimnames = dimnames(lit)
      )
    }
    
    comb <- lit
    comb[!is.na(cor)] <- lit[!is.na(cor)] * cor[!is.na(cor)]
    
    list(
      literature = lit,
      correlation = cor,
      combined = comb
    )
  })
  
  
  # ---- helper: selected matrix and ordering ----
  ct_selected_matrix <- reactive({
    
    req(input$run_crosstalk > 0)
    
    mats <- ct_mats_all()
    
    type <- isolate(input$ct_view_type %||% "combined")
    
    mat <- switch(
      type,
      literature = mats$literature,
      correlation = mats$correlation,
      combined = mats$combined,
      mats$combined
    )
    
    validate(
      need(!is.null(mat), "Cross-talk matrix not available.")
    )
    
    mat_use <- mat
    
    if (!is.null(mat_use)) {
      
      if (
        isolate(input$ct_order_rows) == "Hierarchical clustering" &&
        nrow(mat_use) > 1
      ) {
        d <- stats::dist(mat_use)
        if (all(is.finite(d))) {
          ord <- hclust(d)$order
          mat_use <- mat_use[ord, , drop = FALSE]
        }
      }
      
      if (
        isolate(input$ct_order_cols) == "Hierarchical clustering" &&
        ncol(mat_use) > 1
      ) {
        d <- stats::dist(t(mat_use))
        if (all(is.finite(d))) {
          ord <- hclust(d)$order
          mat_use <- mat_use[, ord, drop = FALSE]
        }
      }
    }
    
    mat_use
  })
  
  
  # ---- reusable heatmap plot function; no numeric labels and no highlight option ----
  make_ct_heatmap_plot <- function(mat) {
    
    df <- as.data.frame(as.table(mat))
    colnames(df) <- c("Ferro", "Sen", "value")
    
    df$Ferro <- factor(df$Ferro, levels = rownames(mat))
    df$Sen   <- factor(df$Sen, levels = colnames(mat))
    
    v <- df$value
    vmax <- max(abs(v[is.finite(v)]), na.rm = TRUE)
    
    if (!is.finite(vmax) || vmax == 0) {
      vmax <- 1
    }
    
    ggplot2::ggplot(
      df,
      ggplot2::aes(x = Sen, y = Ferro, fill = value)
    ) +
      ggplot2::geom_tile(color = "grey80") +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(-vmax, vmax),
        na.value = "grey95",
        name = "Score"
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = "Senescence modules",
        y = "Ferroptosis modules"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          face = "bold"
        ),
        axis.text.y = ggplot2::element_text(face = "bold"),
        axis.title = ggplot2::element_text(face = "bold"),
        panel.border = ggplot2::element_rect(
          colour = "black",
          fill = NA,
          linewidth = 1
        )
      )
  }
  
  
  # ---- Cross-talk heatmap ----
  output$ct_heatmap_matrix <- renderPlot({
    
    req(input$run_crosstalk > 0)
    
    mat <- ct_selected_matrix()
    
    validate(
      need(!is.null(mat), "No cross-talk matrix to display.")
    )
    
    make_ct_heatmap_plot(mat)
  })
  
  
  # ---- Download Cross-talk heatmap PNG ----
  output$ct_heatmap_matrix_png <- downloadHandler(
    
    filename = function() {
      paste0("crosstalk_", input$ct_view_type %||% "combined", ".png")
    },
    
    content = function(file) {
      
      req(input$run_crosstalk > 0)
      
      mat <- ct_selected_matrix()
      p <- make_ct_heatmap_plot(mat)
      
      ggplot2::ggsave(
        file,
        plot = p,
        width = 9,
        height = 7,
        dpi = 300
      )
    }
  )
  
  
  # ---- Download Cross-talk heatmap PDF ----
  output$ct_heatmap_matrix_pdf <- downloadHandler(
    
    filename = function() {
      paste0("crosstalk_", input$ct_view_type %||% "combined", ".pdf")
    },
    
    content = function(file) {
      
      req(input$run_crosstalk > 0)
      
      mat <- ct_selected_matrix()
      p <- make_ct_heatmap_plot(mat)
      
      grDevices::pdf(file, width = 9, height = 7)
      print(p)
      grDevices::dev.off()
    }
  )
  
  
  # ---- Download Cross-talk matrix CSV ----
  output$ct_matrix_csv <- downloadHandler(
    
    filename = function() {
      paste0("crosstalk_", input$ct_view_type %||% "combined", ".csv")
    },
    
    content = function(file) {
      
      req(input$run_crosstalk > 0)
      
      mat <- ct_selected_matrix()
      utils::write.csv(mat, file, row.names = TRUE)
    }
  )
  
  
  # ---- FSI dynamic group selector from metadata ----
  output$fsi_group_col_ui <- renderUI({
    
    meta <- NULL
    
    if (!is.null(aligned())) {
      meta <- tryCatch(aligned()$meta, error = function(e) NULL)
    }
    
    if (is.null(meta) && !is.null(rv$meta)) {
      meta <- rv$meta
    }
    
    if (is.null(meta) || !ncol(as.data.frame(meta))) {
      return(
        selectInput(
          "fsi_group_col",
          "Group samples by metadata column (optional)",
          choices = c("None" = "None"),
          selected = "None"
        )
      )
    }
    
    meta <- as.data.frame(meta)
    choices <- colnames(meta)
    
    preferred <- c(
      "Type", "type",
      "Group", "group",
      "Condition", "condition",
      "Disease", "disease",
      "Treatment", "treatment",
      "Status", "status"
    )
    
    default_choice <- intersect(preferred, choices)
    
    selected_choice <- if (length(default_choice) > 0) {
      default_choice[1]
    } else {
      "None"
    }
    
    selectInput(
      "fsi_group_col",
      "Group samples by metadata column (optional)",
      choices = c("None" = "None", choices),
      selected = selected_choice
    )
  })
  
  
  # ---- FSI data ----
  fsi_data <- reactive({
    
    req(input$run_crosstalk > 0)
    
    sc <- ct_module_scores()
    
    Fm <- sc$ferro
    samples <- sc$samples
    
    validate(
      need(!is.null(Fm), "No module scores available for FSI.")
    )
    
    ferro_protective <- c(
      "Hepatic_Ferroptosis_Core",
      "CoQ10_FSP1_Axis",
      "GCH1_BH4_Axis",
      "NRF2_KEAP1_Antioxidant_Response",
      "Mitochondrial_ROS_and_IronSulfur_Protection"
    )
    
    ferro_proferro <- setdiff(rownames(Fm), ferro_protective)
    
    pro_mat  <- Fm[rownames(Fm) %in% ferro_proferro, , drop = FALSE]
    prot_mat <- Fm[rownames(Fm) %in% ferro_protective, , drop = FALSE]
    
    FSI <- colMeans(pro_mat, na.rm = TRUE) -
      colMeans(prot_mat, na.rm = TRUE)
    
    df <- tibble::tibble(
      Sample = samples,
      FSI = as.numeric(FSI)
    )
    
    meta <- NULL
    
    if (!is.null(aligned())) {
      meta <- tryCatch(aligned()$meta, error = function(e) NULL)
    }
    
    if (is.null(meta) && !is.null(rv$meta)) {
      meta <- rv$meta
    }
    
    if (!is.null(meta) && nrow(meta) > 0) {
      
      meta <- as.data.frame(meta)
      
      # Add Sample column safely.
      if (!"Sample" %in% colnames(meta)) {
        meta$Sample <- rownames(meta)
      }
      
      # If rownames are not useful, use common sample ID columns.
      sample_candidates <- c(
        "Sample", "sample",
        "SampleID", "sample_id",
        "Description", "description",
        "ID", "id"
      )
      
      sample_col <- intersect(sample_candidates, colnames(meta))
      
      if (length(sample_col) > 0) {
        sample_col <- sample_col[1]
        meta$Sample <- as.character(meta[[sample_col]])
      }
      
      df <- dplyr::left_join(df, meta, by = "Sample")
    }
    
    df
  })
  
  
  # ---- FSI color helper ----
  fsi_group_colors <- function(groups) {
    
    groups <- as.character(groups)
    groups <- unique(stats::na.omit(groups))
    
    base_cols <- c(
      "Control" = "#66C2A5",
      "control" = "#66C2A5",
      "CONTROL" = "#66C2A5",
      "Treated" = "#FC8D62",
      "treated" = "#FC8D62",
      "TREATED" = "#FC8D62",
      "Treatment" = "#FC8D62",
      "treatment" = "#FC8D62"
    )
    
    out <- stats::setNames(
      grDevices::hcl.colors(length(groups), palette = "Dark 3"),
      groups
    )
    
    matched <- intersect(names(base_cols), groups)
    out[matched] <- base_cols[matched]
    
    out
  }
  
  
  # ---- FSI density plot helper ----
  fsi_density_plot <- reactive({
    
    req(input$run_crosstalk > 0)
    
    df <- fsi_data()
    req(df)
    
    grp_col <- input$fsi_group_col
    
    if (!is.null(grp_col) &&
        !identical(grp_col, "None") &&
        grp_col %in% names(df)) {
      
      df$Group <- as.factor(df[[grp_col]])
      pal <- fsi_group_colors(levels(df$Group))
      
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = FSI,
          fill = Group,
          color = Group
        )
      ) +
        ggplot2::geom_density(
          alpha = 0.35,
          linewidth = 1.1
        ) +
        ggplot2::scale_fill_manual(values = pal) +
        ggplot2::scale_color_manual(values = pal) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 22,
            face = "bold",
            hjust = 0.5
          ),
          axis.title.x = ggplot2::element_text(
            size = 20,
            face = "bold",
            color = "black",
            margin = ggplot2::margin(t = 10)
          ),
          axis.title.y = ggplot2::element_text(
            size = 20,
            face = "bold",
            color = "black",
            margin = ggplot2::margin(r = 10)
          ),
          axis.text.x = ggplot2::element_text(
            size = 16,
            face = "bold",
            color = "black"
          ),
          axis.text.y = ggplot2::element_text(
            size = 16,
            face = "bold",
            color = "black"
          ),
          legend.title = ggplot2::element_text(
            size = 15,
            face = "bold"
          ),
          legend.text = ggplot2::element_text(
            size = 14,
            face = "bold"
          ),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          x = "FSI",
          y = "Density",
          title = "FSI distribution",
          fill = grp_col,
          color = grp_col
        )
      
    } else {
      
      ggplot2::ggplot(df, ggplot2::aes(x = FSI)) +
        ggplot2::geom_density(
          fill = "grey80",
          color = "black",
          alpha = 0.55,
          linewidth = 1.1
        ) +
        ggplot2::theme_minimal(base_size = 18) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 22,
            face = "bold",
            hjust = 0.5
          ),
          axis.title.x = ggplot2::element_text(
            size = 20,
            face = "bold",
            color = "black",
            margin = ggplot2::margin(t = 10)
          ),
          axis.title.y = ggplot2::element_text(
            size = 20,
            face = "bold",
            color = "black",
            margin = ggplot2::margin(r = 10)
          ),
          axis.text.x = ggplot2::element_text(
            size = 16,
            face = "bold",
            color = "black"
          ),
          axis.text.y = ggplot2::element_text(
            size = 16,
            face = "bold",
            color = "black"
          ),
          panel.grid.minor = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          x = "FSI",
          y = "Density",
          title = "FSI distribution"
        )
    }
  })
  
  
  # ---- FSI by group plot helper ----
  fsi_group_plot <- reactive({
    
    req(input$run_crosstalk > 0)
    req(input$fsi_group_col)
    
    df <- fsi_data()
    req(df)
    
    grp_col <- input$fsi_group_col
    
    if (identical(grp_col, "None")) {
      return(
        ggplot2::ggplot() +
          ggplot2::theme_void() +
          ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "Select a grouping column in metadata.",
            size = 6,
            fontface = "bold"
          ) +
          ggplot2::xlim(0, 1) +
          ggplot2::ylim(0, 1)
      )
    }
    
    validate(
      need(
        grp_col %in% names(df),
        "Selected grouping column not found in metadata."
      )
    )
    
    df$Group <- as.factor(df[[grp_col]])
    
    validate(
      need(
        length(unique(stats::na.omit(df$Group))) >= 1,
        "No valid group labels found for selected metadata column."
      )
    )
    
    pal <- fsi_group_colors(levels(df$Group))
    
    gg <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = Group,
        y = FSI,
        fill = Group,
        color = Group
      )
    )
    
    if (isTRUE(input$fsi_show_violin)) {
      gg <- gg +
        ggplot2::geom_violin(
          alpha = 0.35,
          linewidth = 0.9,
          trim = FALSE
        )
    }
    
    gg +
      ggplot2::geom_boxplot(
        width = 0.25,
        outlier.shape = NA,
        fill = "white",
        color = "black",
        linewidth = 0.9
      ) +
      ggplot2::geom_jitter(
        width = 0.12,
        size = 3,
        alpha = 0.9
      ) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::scale_color_manual(values = pal) +
      ggplot2::theme_minimal(base_size = 18) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = 22,
          face = "bold",
          hjust = 0.5
        ),
        axis.title.x = ggplot2::element_text(
          size = 20,
          face = "bold",
          color = "black",
          margin = ggplot2::margin(t = 10)
        ),
        axis.title.y = ggplot2::element_text(
          size = 20,
          face = "bold",
          color = "black",
          margin = ggplot2::margin(r = 10)
        ),
        axis.text.x = ggplot2::element_text(
          size = 16,
          face = "bold",
          color = "black",
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        axis.text.y = ggplot2::element_text(
          size = 16,
          face = "bold",
          color = "black"
        ),
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = grp_col,
        y = "FSI",
        title = paste0("FSI by ", grp_col)
      )
  })
  
  
  # ---- Render FSI density plot ----
  output$fsi_density <- renderPlot({
    fsi_density_plot()
  })
  
  
  # ---- Render FSI by group plot ----
  output$fsi_by_group <- renderPlot({
    fsi_group_plot()
  })
  
  
  # ---- Download FSI CSV ----
  output$fsi_values_csv <- downloadHandler(
    
    filename = function() {
      "FSI_values.csv"
    },
    
    content = function(file) {
      req(input$run_crosstalk > 0)
      readr::write_csv(fsi_data(), file)
    }
  )
  
  
  # ---- Download current FSI plots as one PDF ----
  output$fsi_current_pdf <- downloadHandler(
    
    filename = function() {
      grp <- input$fsi_group_col %||% "group"
      grp <- gsub("[^A-Za-z0-9]+", "_", grp)
      paste0("FerroEnrich_current_FSI_", grp, "_", Sys.Date(), ".pdf")
    },
    
    content = function(file) {
      
      req(input$run_crosstalk > 0)
      
      p1 <- fsi_density_plot()
      p2 <- fsi_group_plot()
      
      grDevices::pdf(
        file = file,
        width = 12,
        height = 6
      )
      
      grid::grid.newpage()
      grid::pushViewport(
        grid::viewport(
          layout = grid::grid.layout(
            nrow = 1,
            ncol = 2
          )
        )
      )
      
      print(
        p1,
        vp = grid::viewport(
          layout.pos.row = 1,
          layout.pos.col = 1
        )
      )
      
      print(
        p2,
        vp = grid::viewport(
          layout.pos.row = 1,
          layout.pos.col = 2
        )
      )
      
      grDevices::dev.off()
    }
  )
  
  # ---- Top cross-talk edges table and downloads ----
  ct_edges_all <- reactive({
    
    req(input$run_crosstalk > 0)
    
    mats <- ct_mats_all()
    
    lit  <- mats$literature
    cor  <- mats$correlation
    comb <- mats$combined
    
    df <- expand.grid(
      Ferro = rownames(lit),
      Sen = colnames(lit),
      stringsAsFactors = FALSE
    )
    
    df$Literature <- as.vector(lit)
    df$Correlation <- as.vector(cor)
    df$Combined <- as.vector(comb)
    df$absCombined <- abs(df$Combined)
    df$Direction <- ifelse(df$Combined >= 0, "Pro-ferroptotic", "Protective")
    
    df
  })
  
  
  ct_edges_top <- reactive({
    
    req(input$run_crosstalk > 0)
    
    df <- ct_edges_all()
    
    N <- as.integer(input$ct_top_n_edges %||% 30)
    
    if (input$ct_edge_direction_filter == "Pro-ferroptotic only") {
      df <- df[df$Direction == "Pro-ferroptotic", , drop = FALSE]
    } else if (input$ct_edge_direction_filter == "Protective only") {
      df <- df[df$Direction == "Protective", , drop = FALSE]
    }
    
    df <- df[order(-df$absCombined), , drop = FALSE]
    
    df[seq_len(min(N, nrow(df))), , drop = FALSE]
  })
  
  
  output$ct_edges_top_tbl <- renderDT({
    
    req(input$run_crosstalk > 0)
    
    DT::datatable(
      ct_edges_top(),
      options = list(scrollX = TRUE, pageLength = 20),
      rownames = FALSE
    )
  })
  
  
  output$ct_edges_top_csv <- downloadHandler(
    
    filename = function() {
      "Crosstalk_top_edges.csv"
    },
    
    content = function(file) {
      req(input$run_crosstalk > 0)
      readr::write_csv(ct_edges_top(), file)
    }
  )
  
  # ---- Module crosstalk network ----
  output$ct_net_small <- renderPlot({
    
    req(input$run_crosstalk > 0)
    req(.has_igraph)
    
    withProgress(
      message = "Rendering cross-talk network...",
      value = 0,
      {
        incProgress(0.20, detail = "Preparing top module-pair edges...")
        
        edges <- ct_edges_all()
        
        K <- as.integer(input$ct_net_top_k %||% 40)
        
        edges <- edges[order(-edges$absCombined), , drop = FALSE]
        edges <- edges[seq_len(min(K, nrow(edges))), , drop = FALSE]
        
        if (!nrow(edges)) {
          plot.new()
          text(0.5, 0.5, "No edges to plot.", cex = 1.2)
          return()
        }
        
        incProgress(0.25, detail = "Building igraph object...")
        
        edf <- data.frame(
          from = edges$Ferro,
          to = edges$Sen,
          signed_weight = edges$Combined,
          layout_weight = abs(edges$Combined),
          stringsAsFactors = FALSE
        )
        
        edf$layout_weight[is.na(edf$layout_weight)] <- 0
        edf$layout_weight[edf$layout_weight <= 0] <- 1e-6
        
        node_names <- unique(c(edf$from, edf$to))
        
        nodes <- data.frame(
          name = node_names,
          type = ifelse(node_names %in% ct_ferro_mods, "Ferro", "Sen"),
          stringsAsFactors = FALSE
        )
        
        g <- igraph::graph_from_data_frame(
          edf,
          directed = FALSE,
          vertices = nodes
        )
        
        igraph::E(g)$weight <- edf$layout_weight
        igraph::E(g)$signed_weight <- edf$signed_weight
        
        incProgress(0.25, detail = "Calculating network layout...")
        
        set.seed(123)
        
        lay <- switch(
          input$ct_net_layout %||% "fr",
          
          fr = igraph::layout_with_fr(
            g,
            weights = igraph::E(g)$weight
          ),
          
          kk = igraph::layout_with_kk(
            g,
            weights = igraph::E(g)$weight
          ),
          
          circle = igraph::layout_in_circle(g),
          
          igraph::layout_with_fr(
            g,
            weights = igraph::E(g)$weight
          )
        )
        
        # Professional muted node colors
        ferro_col <- "#4E79A7"   # muted blue
        sen_col   <- "#59A14F"   # muted green
        
        vcols <- ifelse(
          igraph::V(g)$type == "Ferro",
          ferro_col,
          sen_col
        )
        
        signed_w <- igraph::E(g)$signed_weight
        
        # Professional muted edge colors
        pro_edge_col  <- "#C44E52"   # muted red
        prot_edge_col <- "#4E79A7"   # muted blue
        
        ecols <- ifelse(
          signed_w >= 0,
          pro_edge_col,
          prot_edge_col
        )
        
        max_abs_w <- max(abs(signed_w), na.rm = TRUE)
        max_abs_w <- ifelse(is.finite(max_abs_w) && max_abs_w > 0, max_abs_w, 1)
        
        ewidth <- 1.2 + 3.2 * abs(signed_w) / max_abs_w
        
        incProgress(0.20, detail = "Drawing cross-talk network...")
        
        par(mar = c(1, 1, 4, 1), xpd = TRUE)
        
        plot(
          g,
          layout = lay,
          
          # Node style
          vertex.color = vcols,
          vertex.size = 30,
          vertex.frame.color = "#2F2F2F",
          vertex.frame.width = 1.2,
          
          # Label style
          vertex.label = igraph::V(g)$name,
          vertex.label.cex = 0.72,
          vertex.label.color = "black",
          vertex.label.font = 1,
          vertex.label.family = "sans",
          vertex.label.dist = 0.45,
          vertex.label.degree = -pi / 4,
          
          # Edge style
          edge.color = grDevices::adjustcolor(ecols, alpha.f = 0.82),
          edge.width = ewidth,
          edge.curved = 0.08,
          
          main = "Ferroptosis-Senescence module crosstalk"
        )
        
        legend(
          "bottomleft",
          legend = c(
            "Ferroptosis module",
            "Senescence module",
            "Pro-ferroptotic edge",
            "Protective edge"
          ),
          pch = c(21, 21, NA, NA),
          pt.bg = c(ferro_col, sen_col, NA, NA),
          pt.cex = c(1.2, 1.2, NA, NA),
          col = c("#2F2F2F", "#2F2F2F", pro_edge_col, prot_edge_col),
          lwd = c(NA, NA, 3, 3),
          bty = "n",
          cex = 0.9
        )
        
        incProgress(0.10, detail = "Network completed.")
        
        gc()
      }
    )
  })
  
  
  # ---- Download Cross-talk network PDF ----
  output$ct_net_pdf_export <- downloadHandler(
    
    filename = function() {
      "Crosstalk_network.pdf"
    },
    
    content = function(file) {
      
      req(input$run_crosstalk > 0)
      req(.has_igraph)
      
      edges <- ct_edges_all()
      
      K <- as.integer(input$ct_net_top_k %||% 40)
      
      edges <- edges[order(-edges$absCombined), , drop = FALSE]
      edges <- edges[seq_len(min(K, nrow(edges))), , drop = FALSE]
      
      if (!nrow(edges)) return(NULL)
      
      edf <- data.frame(
        from = edges$Ferro,
        to = edges$Sen,
        signed_weight = edges$Combined,
        layout_weight = abs(edges$Combined),
        stringsAsFactors = FALSE
      )
      
      edf$layout_weight[is.na(edf$layout_weight)] <- 0
      edf$layout_weight[edf$layout_weight <= 0] <- 1e-6
      
      node_names <- unique(c(edf$from, edf$to))
      
      nodes <- data.frame(
        name = node_names,
        type = ifelse(node_names %in% ct_ferro_mods, "Ferro", "Sen"),
        stringsAsFactors = FALSE
      )
      
      g <- igraph::graph_from_data_frame(
        edf,
        directed = FALSE,
        vertices = nodes
      )
      
      igraph::E(g)$weight <- edf$layout_weight
      igraph::E(g)$signed_weight <- edf$signed_weight
      
      set.seed(123)
      
      lay <- switch(
        input$ct_net_layout %||% "fr",
        
        fr = igraph::layout_with_fr(
          g,
          weights = igraph::E(g)$weight
        ),
        
        kk = igraph::layout_with_kk(
          g,
          weights = igraph::E(g)$weight
        ),
        
        circle = igraph::layout_in_circle(g),
        
        igraph::layout_with_fr(
          g,
          weights = igraph::E(g)$weight
        )
      )
      
      # Professional muted node colors
      ferro_col <- "#4E79A7"
      sen_col   <- "#59A14F"
      
      vcols <- ifelse(
        igraph::V(g)$type == "Ferro",
        ferro_col,
        sen_col
      )
      
      signed_w <- igraph::E(g)$signed_weight
      
      # Professional muted edge colors
      pro_edge_col  <- "#C44E52"
      prot_edge_col <- "#4E79A7"
      
      ecols <- ifelse(
        signed_w >= 0,
        pro_edge_col,
        prot_edge_col
      )
      
      max_abs_w <- max(abs(signed_w), na.rm = TRUE)
      max_abs_w <- ifelse(is.finite(max_abs_w) && max_abs_w > 0, max_abs_w, 1)
      
      ewidth <- 1.2 + 3.2 * abs(signed_w) / max_abs_w
      
      grDevices::pdf(file, width = 11, height = 8.5)
      
      par(mar = c(1, 1, 4, 1), xpd = TRUE)
      
      plot(
        g,
        layout = lay,
        
        # Node style
        vertex.color = vcols,
        vertex.size = 30,
        vertex.frame.color = "#2F2F2F",
        vertex.frame.width = 1.2,
        
        # Label style
        vertex.label = igraph::V(g)$name,
        vertex.label.cex = 0.72,
        vertex.label.color = "black",
        vertex.label.font = 1,
        vertex.label.family = "sans",
        vertex.label.dist = 0.45,
        vertex.label.degree = -pi / 4,
        
        # Edge style
        edge.color = grDevices::adjustcolor(ecols, alpha.f = 0.82),
        edge.width = ewidth,
        edge.curved = 0.08,
        
        main = "Ferroptosis-Senescence module crosstalk"
      )
      
      legend(
        "bottomleft",
        legend = c(
          "Ferroptosis module",
          "Senescence module",
          "Pro-ferroptotic edge",
          "Protective edge"
        ),
        pch = c(21, 21, NA, NA),
        pt.bg = c(ferro_col, sen_col, NA, NA),
        pt.cex = c(1.2, 1.2, NA, NA),
        col = c("#2F2F2F", "#2F2F2F", pro_edge_col, prot_edge_col),
        lwd = c(NA, NA, 3, 3),
        bty = "n",
        cex = 0.9
      )
      
      grDevices::dev.off()
    }
  )
 
  ## =================== end SERVER: Cross-talk tab ===================

}

shinyApp(ui, server)
