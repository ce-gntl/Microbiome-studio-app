#####----------------SERVER-------------------#####
# Define server logic
server <- function(input, output, session) {
  
  #bs_themer()
  
  
  #saves the processed data
  transformed_data <- reactiveValues(OTU = NULL, TAX = NULL, META = NULL,  phy = NULL, norm_phy = NULL,
                                     ancom_results = NULL, ancom_results_norm= NULL)
  
  ########---------INPUT FILES-----------------
  observeEvent(input$process, {
    req(input$otu_file)
    req(input$tax_file)
    req(input$meta_file)
    
    # Read uploaded files
    otu_df <- read.delim(
      input$otu_file$datapath,
      row.names = 1,
      check.names = FALSE
    )
    
    tax_df <- read.delim(
      input$tax_file$datapath,
      row.names = 1,
      sep = ",",
      check.names = FALSE
    )
    
    meta_df <- read.csv(
      input$meta_file$datapath,
      row.names = 1,
      sep = ";",
      check.names = FALSE
    )
    
    # Convert to phyloseq components
    OTU  <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)
    TAX  <- tax_table(as.matrix(tax_df))
    META <- sample_data(meta_df)
    
    # Build phyloseq object
    phy <- phyloseq(OTU, TAX, META)
    
    # Store result
    transformed_data$OTU  <- otu_df
    transformed_data$TAX  <- tax_df
    transformed_data$META <- meta_df
    transformed_data$phy  <- phy
    
    # Populate dropdown choices for tax_level and formula
    updateSelectInput(session, "tax_level", choices = setdiff(colnames(transformed_data$TAX), "Kingdom"))
    updateSelectInput(session, "formula", choices = colnames(transformed_data$META))
    updateSelectInput(session, "meta_column", choices = colnames(transformed_data$META))
    updateSelectInput(session, "group_selection_beta", choices = colnames(transformed_data$META))
    updateSelectInput(session, "comp_group", choices = colnames(transformed_data$META))
    #updateSelectInput(session, "group_selection_alpha", choices = colnames(transformed_data$META))
    
    
    meta_df <- transformed_data$META
    meta <- meta_df[sapply(meta_df, is.character)]
    updateSelectInput(session, "group_selection_alpha", choices = colnames(meta))
  })
  
  #######----------- FILE PREVIEW OUTPUTS---------------
  
  # Output for file preview
  output$preview_table <- renderDT({
    req(transformed_data$OTU, transformed_data$TAX, transformed_data$META)
    
    table_map <- list(
      "Meta Table" = transformed_data$META,
      "Otu Table" = transformed_data$OTU,
      "Taxonomy Table" = transformed_data$TAX
    )
    
    selected_table <- input$preview_choice
    datatable(table_map[[selected_table]], caption = paste(selected_table))
  })
  
  
  ########------------FILES (OTU,TAXA,META) DOWNLOAD------------------
  
  ##download buttons for processed files
  output$downloadOTU<- downloadHandler(
    filename = function() { "otu.csv" },
    content = function(file) { file.copy("otu.csv", file) }
  )
  output$downloadTAXA <- downloadHandler(
    filename = function() { "taxonomy.csv" },
    content = function(file) { file.copy("taxonomy.csv", file) }
  )
  output$downloadMETA <- downloadHandler(
    filename = function() { "meta.csv" },
    content = function(file) { file.copy("meta.csv", file) }
  )
  
  
  
  #########----------ANCOM BC PROCESSING--------------
  
  ##ancom bc processing 
  observeEvent(input$analyze, {
    
    req(transformed_data$phy)
    selected_tax_level <- input$tax_level
    selected_group <- input$formula
    
    
    withProgress(
      min=1, 
      max=10,
      {
        setProgress(1, message = "Here we go")
        out <- ancombc(
          data = transformed_data$phy, tax_level = selected_tax_level, formula = selected_group,
          p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, group = NULL,
          tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, n_cl = 1, verbose = TRUE
        )
        
        transformed_data$ancom_results <- out$res
        transformed_data$out <- out 
        
        setProgress(5, message = "Working hard")
        norm <- ancombc(
          data = transformed_data$phy, tax_level = NULL, formula = selected_group,
          p_adj_method = "holm", prv_cut = 0.10, lib_cut = 0, group = NULL,
          tol = 1e-5, max_iter = 100, conserve = TRUE, alpha = 0.05, n_cl = 1, verbose = TRUE
        )
        
        transformed_data$ancom_results_norm <- norm$res
        transformed_data$norm <- norm
        
        setProgress(10, message = "Finished!")
      }
    )
    
    
    ##making the new normalized dataframe with abundances
    samp_frac = norm$samp_frac
    # Replace NA with 0
    samp_frac[is.na(samp_frac)] = 0
    # Add pesudo-count (1) to avoid taking the log of 0
    log_obs_abn = log(norm$feature_table + 1)
    # Adjust the log observed abundances
    log_corr_abn = t(t(log_obs_abn) - samp_frac)
    #print(log_corr_abn)
    
    OTU <- otu_table(as.matrix(log_corr_abn), taxa_are_rows= TRUE)
    #print(OTU)
    TAX <- tax_table(as.matrix(transformed_data$TAX))
    #print(TAX)
    META <- sample_data(transformed_data$META)
    norm_phy <- phyloseq(OTU, TAX, META)
    transformed_data$norm_phy <- norm_phy
    
  })
  
  #######-----------ANCOM BC RESULTS OUTPUT---------------------
  
  # Output for ANCOMBC results
  output$ancombc_table <- renderDT({
    req(transformed_data$ancom_results)
    
    
    res <- transformed_data$ancom_results
    
    table_map <- list(
      "Log Fold Changes" = res$lfc,
      "Standard Errors" = res$se,
      "Test Statistics" = res$W,
      "P-values" = res$p_val,
      "Adjusted P-values" = res$q,
      "Differentially Abundant Taxa" = res$diff_abn,
      "Bias-corrected Abundances" = {
        req(transformed_data$out)
        samp_frac <- transformed_data$out$samp_frac
        samp_frac[is.na(samp_frac)] <- 0
        log_obs_abn <- log(transformed_data$out$feature_table + 1)
        log_corr_abn <- t(t(log_obs_abn) - samp_frac)
        log_corr_abn[, 1:6]
      }
    )
    
    selected_table <- input$result_choice
    datatable(table_map[[selected_table]], caption = paste(selected_table, "from the Primary Result"))
  })
  
  
  ##########-----------DATAFRAME CREATION FOR LFC DATASET AND PLOT--------------
  
  # Compute significant taxa and prepare LFC DataFrame
  df_fig_group <- reactive({
    req(transformed_data$ancom_results)
    col_name = c("Taxon", "Intercept", "Group")
    
    res <- transformed_data$ancom_results
    
    tab_lfc = res$lfc
    colnames(tab_lfc) = col_name
    tab_lfc %>% 
      datatable() %>%
      formatRound(col_name[-1], digits = 5)
    
    tab_se <- res$se
    colnames(tab_se) = col_name
    tab_se %>% 
      datatable() %>%
      formatRound(col_name[-1], digits = 5)
    
    tab_p <- res$q
    colnames(tab_p) = col_name
    tab_p %>% 
      datatable() %>%
      formatRound(col_name[-1], digits = 5)
    
    ###slects only significant taxa (with p_val<0,05)
    sig_taxa <- tab_p %>%
      filter(tab_p[,3]< 0.05) %>%
      .$Taxon
    ###filters lfc of only significant taxa (as selcted in previous step)
    df_group = tab_lfc %>%
      dplyr::select(Taxon, Group) %>%
      filter(Taxon %in% sig_taxa)
    
    ##filters SE of selected taxa
    df_se = tab_se %>%
      dplyr::select(Taxon, Group) %>%
      filter(Taxon %in% sig_taxa)
    colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "(SE)")
    
    ###filters p-values of selected taxa
    df_p = tab_p %>%
      dplyr::select(Taxon, Group) %>%
      filter(Taxon %in% sig_taxa)
    colnames(df_p)[-1] <- c("p_value")
    
    ######Arranges df for plotting (by decreasing fold change value)
    df_fig_group = df_group %>% 
      transmute(Taxon, 
                `Group comparison` = round(Group, 5)) %>%
      pivot_longer(cols =  `Group comparison`, 
                   names_to = "group", values_to = "value") %>%
      arrange(value)
    
    #adds one extra column to dataset, assign labels to positive and negative fold change
    direct = ifelse(df_fig_group$value > 0, "Positive LFC", "Negative LFC")
    df_fig_group['direct'] <- direct
    #df_fig_group$direct = factor(df_fig_group$direct)
    
    ###adds SE to the fold change dataframe
    df_fig_group = df_fig_group %>% 
      dplyr::inner_join(df_se, by = "Taxon")
    
    ###adds p_val to the fold change dataframe
    df_fig_group = df_fig_group %>% 
      dplyr::inner_join(df_p, by = "Taxon")
    df_fig_group$group <- NULL
    
    return(df_fig_group)
  })
  
  
  ######---------DATAFRAME LFC----------------
  
  # Render DataTable for LFC values
  output$lfc_table <- renderDT({
    req(df_fig_group())
    datatable(df_fig_group(), caption = "Log Fold Change Data")
  })
  
  ########----------FOLD CHANGE PLOT---------------
  
  # Render Fold Change Plot
  output$FoldChange <- renderPlot({
    req(df_fig_group())
    
    ggplot(data = df_fig_group(), aes(x = reorder(Taxon, value), y = value, fill = direct, color = direct)) + 
      geom_bar(stat = "identity", width = 0.7, position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin = value - `Group(SE)`, ymax = value + `Group(SE)`), width = 0.2,
                    position = position_dodge(0.05), color = "gray") + 
      labs(x = NULL, y = "Log fold change", title = "Log fold changes") + 
      scale_fill_discrete(name = NULL) +
      scale_color_discrete(name = NULL) +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1))
  })
  
  ######-----------VENN DIAGRAM---------------
  
  # Dynamic UI for group selection in venn diagram
  output$group_select_ui <- renderUI({
    req(input$meta_column, transformed_data$META)
    
    unique_groups <- unique(transformed_data$META[[input$meta_column]])
    
    if (length(unique_groups) < 2) {
      return(NULL)
    }
    
    selectInput("selected_groups", "Select up to ten Groups", choices = unique_groups, multiple = TRUE, selected = unique_groups[1:2])
  })
  
  # Generate Venn Diagram Data
  venn_data <- eventReactive(input$generate_venn, {
    req(transformed_data$phy, input$meta_column, input$selected_groups)
    
    # Extract data from phyloseq object
    bac_meta_df <- psmelt(transformed_data$phy)
    
    # Validate group selection
    if (length(input$selected_groups) > 10) {
      showNotification("Please select maximum ten groups.", type = "error")
      return(NULL)
    }
    
    # Create Venn data
    venn_data <- list()
    
    for (i in seq_along(input$selected_groups)) {
      group_name <- input$selected_groups[i]
      venn_data[[group_name]] <- unique(bac_meta_df$OTU[
        bac_meta_df[[input$meta_column]] == group_name & bac_meta_df$Abundance > 0
      ])
    }
    
    return(venn_data)
  })
  
  
  #######-------------VENN DIAGRAM OUTPUT---------------
  # Render Venn Diagram
  output$venn_diagram <- renderPlot({
    req(venn_data())
    ggvenn(venn_data(), stroke_size = 0.5, fill_color = c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31)))
  })
  
  
  
  ##########--------------COMPOSITIONAL ANALYSIS-----------------------------
  
  output$selected_plot <- renderPlot({
    req(transformed_data$phy, input$comp_group, input$tax_level, input$detection, input$prevalence)
    
    phy<- transformed_data$phy
    norm_phy <- transformed_data$norm_phy
    
    #print(colnames(as.data.frame(tax_table(phy))))
    # Check if input$tax_level exists in phy's tax_table
    if (!(input$tax_level %in% colnames(as.data.frame(tax_table(phy))))) {
      stop(paste("Error: Column", input$tax_level, "not found in tax_table."))
    }
    
    mycols = c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))
    
    Group <- toString(input$comp_group)
    # # Transform to compositional abundances
    # pseq.rel <- microbiome::transform(phy, "compositional")
    
    
    # Filter by taxonomic level, detection, and prevalence
    # pseq.filtered <- aggregate_rare(phy, level = input$tax_level, 
    #                                 detection = input$detection / 100, 
    #                                 prevalence = input$prevalence / 100)
    
    
    #####-----core composition heatmap------------------------------
    
    if (input$plot_choice == "Core Composition Heatmap") {
      pseq.rel <- microbiome::transform(phy, "compositional")
      prevalences <- seq(.05, 1, .05)
      detections <- round(10^seq(log10(5e-3), log10(.2), length = 10), 3)
      p <- plot_core(pseq.rel, plot.type = "heatmap",
                     prevalences = prevalences, detections = detections, min.prevalence = 0.5) +
        xlab("Detection Threshold (Relative Abundance)") +
        theme(axis.text.x = element_text(size = 9))
      return(p)
    }
    
    #######----------------Bar plot group differences------------------------anc
    
    else if (input$plot_choice == "Order Bar Plot") {
      
      minTotRelAbun = 1e-2
      x = taxa_sums(phy)
      keepTaxa = (x / sum(x)) > minTotRelAbun
      prunedSet = prune_taxa(keepTaxa, phy)
      
      grp_abund <- get_group_abundances(prunedSet, 
                                        level = input$tax_level, 
                                        group=Group,
                                        transform = "compositional")
      
      
      mean.plot.ord <- grp_abund %>% # input data
        ggplot(aes(x= reorder(OTUID, mean_abundance), # reorder based on mean abundance
                   y= mean_abundance,
                   fill=!!sym(Group))) + # x and y axis 
        geom_bar(stat = "identity",
                 position=position_dodge()) + 
        #scale_fill_manual(values=mycols) + # manually specify colors
        theme_bw() + # add a widely used ggplot2 theme
        ylab("Mean Relative Abundance") + # label y axis
        xlab(input$tax_level) + # label x axis
        coord_flip() # rotate plot 
      
      return(mean.plot.ord)
    }
    
    ####--------------- Pie Charts----------------
    
    else if (input$plot_choice == "Pie Charts") {
      
      req(input$comp_group, input$tax_level)
      
      
      ## Agglomerate at selected taxonomic rank
      phy_glom <- tax_glom(phy, taxrank = input$tax_level)
      
      ## Transform to relative abundance (%)
      phy_rel <- transform_sample_counts(phy_glom,function(x) x / sum(x) * 100)
      
      phy_melt <- psmelt(phy_rel)
      
      phy_avg <- phy_melt %>%
        filter(!is.na(Abundance)) %>%
        group_by(
          Group = .data[[input$comp_group]],
          Taxon = .data[[input$tax_level]]
        ) %>%
        summarise(
          Abundance = mean(Abundance, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        filter(Abundance > 0) %>%
        group_by(Group) %>%
        slice_max(Abundance, n = 10, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(Percentage = round(Abundance, 1)) %>%
        arrange(Group, desc(Abundance)) %>%
        group_by(Group) %>%
        mutate(
          ypos = cumsum(Abundance) - Abundance / 2
        )
      
      
      ## Plot
      pie_plot <- ggplot(phy_avg, aes(x = "", y = Abundance, fill = Taxon)) +
        geom_col(width = 1) +
        facet_wrap(~ Group) +
        coord_polar("y") +
        theme_void() +
        theme(
          legend.position = "bottom"
        ) +
        labs(fill = input$tax_level) +
        scale_fill_brewer(palette = "Set3") +
        geom_label_repel(
          aes(
            x = 1.6,
            label = paste0(Percentage, "%")
          ),
          size = 3,
          show.legend = FALSE,
          color = "black",
          position = position_stack(vjust = 0.5)
        )
      
      return(pie_plot)
      
    }
    ########---------------Violin plots----------------------------#######
    
    else if (input$plot_choice == "Relative Abundance Plot by Patient") {
      Group <- "Group"
      
      phylumabundance <- phy %>%
        tax_glom(taxrank = "Family") %>%                     # Agglomerate at family or genus level
        transform_sample_counts(function(x) {x / sum(x)}) %>% # Transform to rel. abundance
        psmelt() %>%                                         # Melt to long format
        arrange(Family) 
      
      # Select the top 10 most abundant families
      top_families <- phylumabundance %>%
        group_by(Family) %>% ##change for family or genus
        summarise(mean_abundance = mean(Abundance)) %>%
        arrange(desc(mean_abundance)) %>%
        slice_head(n = 10) %>%
        pull(Family) ##change for family or genus
      
      # Filter for only the top 10 families
      phylumabundance <- phylumabundance %>%
        filter(Family %in% top_families) ##change for family or genus
      
      label_trueminus <- function(x){
        ifelse(sign(x) == -1, paste0("\u2212", abs(x)), x)
      }
      
      # Plot
      library(scales)
      ggplot(phylumabundance, aes(x = Group, y = Abundance, fill = Group)) +
        geom_violin(trim = FALSE, alpha = 0.7) + 
        geom_jitter(width = 0.2, alpha = 0.2, color = "black") +  # Add points for visibility
        facet_wrap(~Family, scales = "free", ncol = 4) +  # Separate plots for each family ##change for family or genus
        theme_minimal() +
        labs(title = "Top 10 Most Abundant Families Across Groups", ##change for family or genus
             x = "Group",
             y = "Relative Abundance") +
        #scale_y_continuous(labels = label_comma())+
        scale_y_continuous(labels = label_trueminus)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size= 15),
              strip.text.x = element_text(size = 15))
    }
    
    
    ######-----------------Relative Abundance Plot by Patient------------------
    
    else if (input$plot_choice == "Relative Abundance Plot by Patient") {
      pseq.filtered <- aggregate_rare(phy, level = input$tax_level, 
                                      detection = input$detection / 100, 
                                      prevalence = input$prevalence / 100)
      pseq.rel <- microbiome::transform(pseq.filtered, "compositional")
      
      p.famrel <- plot_composition(pseq.rel,
                                   sample.sort = NULL,
                                   otu.sort = NULL,
                                   x.label = "empo_3",
                                   plot.type = "barplot",
                                   verbose = FALSE) +
        theme_minimal() +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "Samples",
             y = "Relative abundance",
             title = "Relative abundance data",
             subtitle = "",
             caption = "") +
        scale_fill_brewer(input$tax_level, palette = "Set3") +
        #theme(axis.text.x = element_text(angle=180, hjust=1),
        #     legend.text = element_text(face = "italic"))
        #Removes sample names and ticks
        
        theme(axis.text.x=element_blank(), 
              axis.ticks.x=element_blank()) +
        #Adjusts size of subtitle, caption, legend text and legend title
        theme(plot.subtitle=element_text(size=8), 
              plot.caption=element_text(size=8), 
              legend.text=element_text(size=8),
              legend.title =element_text(size=10))
      
      return(p.famrel)
    }
    
    ######-------------------Relative Abundance Divided by Groups----------------
    
    else if (input$plot_choice == "Relative Abundance Divided by Groups") {
      pseq.filtered <- aggregate_rare(phy, level = input$tax_level, 
                                      detection = input$detection / 100, 
                                      prevalence = input$prevalence / 100)
      # Transform counts to relative abundance
      ps.rel <- transform_sample_counts(pseq.filtered, function(x) x / sum(x) * 100)
      #ps.rel <- microbiome::transform(pseq.filtered, "compositional")
      
      # Agglomerate taxa
      glom <- tax_glom(ps.rel, taxrank = input$tax_level, NArm = FALSE)
      ps.melt <- psmelt(glom)
      
      # Debug: Print available columns
      #print(colnames(ps.melt))  
      #print(paste("Grouping by:", input$comp_group))  
      
      # Dynamically reference the grouping column
      group_col <- sym(input$comp_group)
      
      # Check if the column exists
      if (!(input$comp_group %in% colnames(ps.melt))) {
        stop(paste("Error: Column", input$comp_group, "not found in data."))
      }
      
      # Convert tax_level to character
      tax_col <- sym(input$tax_level)
      ps.melt <- ps.melt %>%
        mutate(!!tax_col := as.character(!!tax_col))
      
      # Compute median abundance per group
      ps.melt <- ps.melt %>%
        group_by(!!group_col, !!tax_col) %>%
        mutate(median = median(Abundance, na.rm = TRUE))
      
      # Select taxa with median abundance > 2.5%
      keep <- unique(ps.melt[[input$tax_level]][ps.melt$median > 2.5])
      ps.melt[[input$tax_level]][!(ps.melt[[input$tax_level]] %in% keep)] <- "< 2.5%"
      
      # Summarize abundance by Sample, group, and tax level
      ps.melt_sum <- ps.melt %>%
        group_by(Sample, !!group_col, !!tax_col) %>%
        summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop')
      
      # Plot
      p <- ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = !!tax_col)) + 
        geom_bar(stat = "identity") + 
        facet_wrap(vars(!!group_col), scales = "free_x", nrow = 1) +
        theme_minimal() +
        guides(fill = guide_legend(ncol = 1)) +
        labs(x = "Samples",
             y = "Relative abundance",
             title = "Relative Abundance Divided by Groups") +
        scale_fill_brewer(palette = "Set3") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.text = element_text(face = "italic", size = 8),
              plot.subtitle = element_text(size = 8), 
              plot.caption = element_text(size = 8), 
              legend.title = element_text(size = 10))
      
      return(p)
      
    }
  })
  
  
  
  ########----------ALPHA DIVERSITY PLOT---------------------
  
  observeEvent(input$alpha_process, {
    req(transformed_data$norm_phy, input$formula, input$group_selection_alpha)
    
    group <-input$group_selection_alpha
    set.seed(4235421)
    # Define color palette
    mycols = c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))
    
    # Generate Alpha Diversity Plots
    obs_alpha_plot <- plot_diversity_stats(transformed_data$norm_phy, group = group, 
                                           index = "observed", label.format="p.format", 
                                           group.colors = mycols, stats = TRUE)
    
    chao1_alpha_plot <- plot_diversity_stats(transformed_data$norm_phy, group = group, 
                                             index = "chao1", label.format="p.format", 
                                             group.colors = mycols, stats = TRUE)
    
    shan_alpha_plot <- plot_diversity_stats(transformed_data$norm_phy, group = group, 
                                            index = "diversity_shannon", label.format="p.format", 
                                            group.colors = mycols, stats = TRUE)
    
    invsimp_alpha_plot <- plot_diversity_stats(transformed_data$norm_phy, group = group, 
                                               index = "diversity_inverse_simpson", label.format="p.format", 
                                               group.colors = mycols, stats = TRUE)
    
    # Arrange plots in a grid
    alphadiv_wp <- ggarrange(obs_alpha_plot, chao1_alpha_plot, 
                             shan_alpha_plot, invsimp_alpha_plot, 
                             ncol = 4, nrow = 1)
    
    # Store plots in reactiveValues
    transformed_data$alpha_plot <- alphadiv_wp
    
  })
  
  
  #############---------------alpha statistics------------------------------
  observeEvent(input$alpha_process,{
    req(transformed_data$norm_phy, input$group_selection_alpha)
    
    group <- input$group_selection_alpha
    ps <- transformed_data$norm_phy
    
    # Construct the data
    d <- meta(ps)
    tab <-microbiome::alpha(ps, index = "all")
    indexes <- list("chao1", "observed","diversity_shannon","diversity_inverse_simpson")
    
    # Create an empty dataframe to store all results
    final_results <- data.frame(
      Index = character(),
      Group1 = character(),
      Group2 = character(),
      P_Value = numeric(),
      Adjusted_P_Value = numeric(),
      stringsAsFactors = FALSE
    )
    
    for(index in indexes) {
      
      d$diversity <- tab[[index]]
      
      # Get unique groups in the grouping column
      unique_groups <- unique(d[[group]])
      
      # Create a dataframe to store results for this index
      index_results <- data.frame(
        Index = character(),
        Group1 = character(),
        Group2 = character(),
        P_Value = numeric(),
        Adjusted_P_Value = numeric(),
        stringsAsFactors = FALSE
      )
      
      # Perform KS test for each pair of unique values
      for (i in 1:(length(unique_groups) - 1)) {
        for (j in (i + 1):length(unique_groups)) {
          group1 <- unique_groups[i]
          group2 <- unique_groups[j]
          
          spl1 <- d$diversity[d[[group]] == group1]
          spl2 <- d$diversity[d[[group]] == group2]
          
          ks_result <- ks.test(spl1, spl2)
          
          # Store results in dataframe
          index_results <- rbind(index_results, data.frame(
            Index = index,
            Group1 = group1,
            Group2 = group2,
            P_Value = ks_result$p.value,
            Adjusted_P_Value = NA
          ))
        }
      }
      
      # Adjust the p-values for multiple comparisons
      if (nrow(index_results) > 0) {
        index_results$Adjusted_P_Value <- p.adjust(index_results$P_Value)
      }
      final_results <- rbind(final_results, index_results)
    }
    transformed_data$alpha_statistics <- final_results
  }
  
  )
  output$alpha_statistics <- renderDT({
    req(transformed_data$alpha_statistics)
    transformed_data$alpha_statistics
  })
  
  
  #render alpha diversity data table
  output$alpha_index <- renderDT({
    req(transformed_data$norm_phy)
    tab <-microbiome::alpha(transformed_data$norm_phy, index = "all")
  })
  
  
  # Render Alpha Diversity Plot
  output$alpha_plots <- renderPlot({
    req(transformed_data$alpha_plot)
    transformed_data$alpha_plot
  })
  
  
  ############-------------------BETA DIVERSITY----------------------------------------
  
  ###---pcoa plot----
  
  output$pcoa_plot <- renderPlot({
    req(transformed_data$norm_phy, input$beta_index, input$group_selection_beta)
    
    phy <- transformed_data$norm_phy
    index <- toString(input$beta_index) 
    group <- input$group_selection_beta
    
    set.seed(4235421)
    #PCoA with jaccard index
    ord.pcoa <- ordinate(phy, "PCoA", index)
    ord_PCoA = plot_ordination(phy, ord.pcoa, color = group) +
      scale_fill_distiller(palette = "Set2")+
      geom_point(size = 2)+ stat_ellipse()+
      theme_classic()
    
    
    return(ord_PCoA)
  })
  
  ####---permanova-----
  
  output$permanova_out <- renderPrint({
    req(transformed_data$norm_phy, input$beta_index, input$group_selection_beta, input$beta_index)
    
    set.seed(4235421)
    
    phy <- transformed_data$norm_phy
    #index <- toString(input$beta_index) 
    group <- input$group_selection_beta
    index <- toString(input$beta_index)
    
    # Generate distance matrix using Bray-Curtis method
    dist_matrix <- phyloseq::distance(phy, method = index) 
    
    # Run ADONIS test with dynamically selected grouping variable
    perm <- vegan::adonis2(dist_matrix ~ phyloseq::sample_data(phy)[[group]])
    
    # Return formatted output for DataTable
    return(perm)
  })
  
  
  ###---dispersion-----
  
  output$dispr <- renderPrint({
    req(transformed_data$norm_phy, input$beta_index, input$group_selection_beta, input$beta_index)
    
    set.seed(4235421)
    
    phy <- transformed_data$norm_phy
    #index <- toString(input$beta_index) 
    group <- input$group_selection_beta
    index <- toString(input$beta_index)
    
    # Generate distance matrix using Bray-Curtis method
    dist_matrix <- phyloseq::distance(phy, method = index) 
    dispr <- vegan::betadisper(dist_matrix, phyloseq::sample_data(phy)[[group]])
    transformed_data$dispr <- dispr
    
    return(dispr)
  })
  
  output$disp_plot <- renderPlot({
    req(transformed_data$dispr, input$disp_plot_sel)
    
    dispr <- transformed_data$dispr
    
    if (input$disp_plot_sel == "Ordination Centroids and Dispersion Labeled") {
      disp_plot <- plot(dispr, main = "Ordination Centroids and Dispersion Labeled", sub = "")
      
      return(disp_plot)
    }
    
    else {
      box <- boxplot(dispr, main = "", xlab = "")
      
      return(box)
    }
    
  })
  
  output$permutest <- renderPrint({
    req(transformed_data$dispr)
    
    x <- permutest(transformed_data$dispr)
    return(x)
  })
  
}
