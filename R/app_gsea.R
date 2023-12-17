#' Shiny APP GSEA
#' @description Shiny APP GSEA
#' @param gmt.largelist dataframe or gmt file from msigdb
#' @param genelist list, a name list, gene list sorted by log2FC 
#' 
#' @return no
#' @import shiny
#'
#' @author Jiang
syGSEA <- function(gmt.largelist, genelist) {
  source.name <- names(gmt.largelist)
  shinyApp(
    ui = fluidPage(
      tags$head(
        tags$title("GSEA"),
        tags$style(
          HTML(
            ".panel-body {
              padding-top: 8px;
            }
            .btn-primary {
              background-color: #3498db;
              border-color: #3498db;
            }
            .btn-primary:hover {
              background-color: #2980b9;
              border-color: #2980b9;
            }
            .panel-primary > .panel-heading {
              color: #fff;
              background-color: #3498db;
              border-color: #3498db;
            }
            .selectize-dropdown {
              z-index: 9999;
            }
            .notification-center {
            position: fixed;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            background-color: rgba(128, 128, 128, 0.1);
            color: rgba(70,130,180,1);
            font-weight: bold;
            font-size: 24px;
            position: fixed;
            padding: 8px;
            border-radius: 4px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
            z-index: 9999;
            }"
          )
        )
      ),
      # tags$div(h4("Gene Set Enrichment Analysis(GSEA)")),
      tags$div(HTML("<p style='font-size: 18px; font-family: Microsoft Yahei; font-weight:bold; color: #3498db;'> Gene Set Enrichment Analysis(GSEA) </p>")),
      fluidRow(
        column(
          width = 6,
          tags$div(
            class = "panel panel-primary",
            tags$div(
              class = "panel-body",
              selectInput("source_set", 
                          label = "通路数据集来源选择",
                          choices = source.name, 
                          selected = NULL,
                          width = '100%'),
              selectInput("pathway_name", 
                          label = "通路名称：选择或输入关键字查找",
                          choices = NULL, 
                          width = '100%'),
              actionButton("plot_btn", "画图", class = "btn-primary"),
              downloadButton("download_plot", "下载", class = "btn-primary"),
              actionButton("submit_btn", "结束", class = "btn-primary")
            )
          ),
          br(),
          tags$div(
            class = "panel panel-primary",
            tags$div(
              class = "panel-body",
              uiOutput('pathway_title'),
              textOutput('pathway_genelist')
            )
          ),
          fileInput("upload_file", "上传表格文件（XLSX格式）", accept = ".xlsx")
        ),
        column(
          width = 6,
          plotOutput("plot_gsea")
        )
      ),
      fluidRow(
        column(
          width = 6,
          
        )
      )
    ),
    server = function(input, output, session) {
      selected_pathway <- reactiveVal(NULL)  # 创建一个响应式对象来存储选择的值
      gmt_data <- reactiveValues(gmt = NULL, gmt_taget = NULL)  # 创建一个reactiveValues对象来存储gmt数据
      plot_requested <- reactiveVal(FALSE)
      # Server logic for handling file upload
      uploaded_data <- reactiveVal(NULL)
      xlsx_data <- reactiveValues(genelist = genelist)
      pic <- reactiveValues(gsea = NULL)
      
      
      observeEvent(input$source_set, {
        # 获取第一个选择的值
        selected_option <- input$source_set
        gmt_data$gmt <-  gmt.largelist[[selected_option]]
        source.path.name <- sort(unique(gmt_data$gmt$term))
        
        # 根据第一个选择更新第二个选择的内容
        updateSelectInput(session, "pathway_name",
                          choices = source.path.name)
      })
      
      
      output$pathway_title <- renderUI({
        if (!is.null(gmt_data$gmt) & !is.null(input$pathway_name)) {
          HTML(paste0("<div style='font-size: 18px; font-family: Arial; font-weight:bold; color: #FF0000;'>", 
                      input$pathway_name, " 通路基因列表：",
                      "</div>")
          )
        }
      })
      
      
      output$pathway_genelist <- renderText({
        if (!is.null(gmt_data$gmt) & !is.null(input$pathway_name)) {
          gene <- gmt_data$gmt[gmt_data$gmt$term == input$pathway_name, 2]
          c(sapply(gene[1:(length(gene)-1)], function(x) paste0(x, ",")), 
            tail(gene, 1))
        }
      })
      
      
      observeEvent(input$plot_btn, {
        tryre <- tryCatch({
          if (length(xlsx_data$genelist) == 1) {
            showModal(modalDialog(
              title = "提示",
              HTML("<p>由于启动时未提供genelist，因此需要用户上传genelist数据。两种解决方法:
            <br>1. 请按要求上传数据
            <br>2. 或关闭程序在下次启动时设置genelist参数</p>"),
              easyClose = TRUE,
              footer = modalButton("关闭")
            ))
          } else {
            id <- showNotification('通知', 
                                   tags$div(
                                     "正在计算中，请勿关闭网页...",
                                     class = "notification-center"
                                   ),
                                   duration = NULL, closeButton = FALSE)
            on.exit(removeNotification(id), add = TRUE)
            taget_gmtdf <- gmt_data$gmt[gmt_data$gmt$term == input$pathway_name, ]
            gsea <- DEG_runGSEA(genelist=xlsx_data$genelist, gmt_set=taget_gmtdf, pic.save=F)
            gmt_data$gmt_taget <- taget_gmtdf
            pic$gsea <- DEGp_GSEA(gsea, num = 1)
            pic$data <- gsea
            # Sys.sleep(10)
            output$plot_gsea <- renderPlot(width = 600, height=514, res = 100,
                                           pic$gsea %>% print() )
          }
        }, error = function(err) {
          stop(cat(conditionMessage(err), "\n"))
        })
      })
      
      
      observeEvent(input$upload_file, {
        uploaded_data(NULL)
        xlsx_data$genelist <- NA
        req(input$upload_file) # Require the file to be uploaded
        uploaded <- readxl::read_xlsx(input$upload_file$datapath, sheet = 1)
        
        # Check if uploaded file has two columns: Gene and log2FC
        if (ncol(uploaded) == 2 && all(c("Gene", "log2FC") == colnames(uploaded))) {
          uploaded_data(uploaded)
        } else {
          uploaded_data(NULL)
          showModal(modalDialog(
            title = "文件检验失败",
            "请确保上传的文件仅包含两列，第一列为Gene，第二列为log2FC。",
            easyClose = TRUE,
            footer = modalButton("关闭")
          ))
        }
        
        if (!is.null(uploaded_data())) {
          genelist <- uploaded_data()$log2FC # Assigning the Gene column to genelist
          names(genelist) <- uploaded_data()$Gene
          genelist <- sort(genelist, decreasing = T)
          xlsx_data$genelist <- genelist
        } else {
          updateFileInput(session, "upload_file", "上传表格文件（XLSX格式）", accept = ".xlsx")
        }
      })
      
      output$download_plot <- downloadHandler(
        filename = function() {
          if (is.null(pic$gsea)) { 
            "please plot picture and then click download button.txt"
          } else {
            paste0('gsea_plot_', input$pathway_name, 
                   format(Sys.time(), "-%Y-%m-%d-%H%:%M:%S"), '.pdf')
          }
        },
        content = function(file) {
          if (is.null(pic$gsea)) { 
            writeLines("错误：请先按要求画图后再点击下载按键", con = file) 
          } else {
            pdf(file, width = 7, height = 6)
            print(pic$gsea)
            dev.off()
          }
        }
      )
      
      # 关闭应用时的动作
      observeEvent(input$submit_btn, {
        selected_pathway(input$pathway_name)  # 将选择的值存储到响应式对象中
        # 关闭Shiny应用程序
        session$reload()  # 重新加载应用程序，模拟关闭应用的效果
      })
      # 当应用程序重新加载时，执行一些操作
      observe({
        if (!is.null(selected_pathway())) {
          cat("Selected Pathway is:", selected_pathway())
          assign("sy_pathway_name", selected_pathway(), envir = .GlobalEnv)
          assign("sy_gmt", gmt_data$gmt, envir = .GlobalEnv)
          assign("sy_gmt_taget", gmt_data$gmt_taget, envir = .GlobalEnv)
          assign("pic_gsea", pic$gsea, envir = .GlobalEnv)
          stopApp()  # 停止Shiny应用程序
        }
      })
      # Function to execute when the session ends
      onSessionEnded(function() {
        # Clean-up operations or actions when the session ends
        stopApp()  # Stop the Shiny app explicitly
      })
      
    }
  )
}


#' Run Shiny APP GSEA
#' @description Run Shiny APP GSEA
#' @param genelist list, a name list, gene list sorted by log2FC 
#' @return no
#' @export
#' @import shiny
#'
#' @author Jiang
runAPP_GSEA <- function(genelist = NA) {
  library(LZ)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(shiny)
  tmpenv <- new.env(parent = emptyenv())
  data(gmt.largelist.23.12.Hs.symbols, package = "LZ", envir = tmpenv)
  gmt.largelist <- tmpenv$gmt.largelist.23.12.Hs.symbols
  rm(tmpenv)
  syGSEA(gmt.largelist = gmt.largelist, genelist = genelist)
}

