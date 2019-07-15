library(shiny)
library(ggplot2)

options(shiny.maxRequestSize=30*1024^2,shiny.trace=F,shiny.fullstacktrace=F,shiny.testmode=F)

function(input,output,session) {

	##  initialise reactive values here
	names <- reactiveValues(
		bed = NULL,
		plot = NULL,
		currentRatio = NULL,
		ratio = data.frame(name.rep=character(),name.nonRep=character(),name=character(),stringsAsFactors=F)
	)
	DFs <- reactiveValues(
		bed = NULL,
		rep = list(),
		nonRep = list(),
		ratio = NULL,
		ratios2add = data.frame(
			name.rep=character(),
			name.nonRep=character()
		),
		ratios = NULL,
		guide = data.frame(
			order=integer(),
			name.rep=character(),
			name.nonRep=character(),
			raw=logical(),
			smooth=logical(),
			color=character(),
			stringsAsFactors=F
		),
		stats = NULL,
		lines = NULL,
		circles = NULL,
		rectangles = NULL,
		pointers = NULL,
		Lines = NULL,
		Circles = NULL,
		Rectangles = NULL,
		Pointers = NULL
	)
	plots <- reactiveValues(coveragePlot = NULL, ratioPlot = NULL, genomePlot = NULL, statsPlot = NULL)
	toggles <- reactiveValues(
		bedCtrlsIn = F,
		bedEdited = F,
		coverage = F,
		isExampleCoverage = F,
		ratioCtrlsIn = F,
		isExampleRatio = F,
		plotCtrlsIn = F,
		isExamplePlot = F,
		isExampleStats = F,
		coverageRegion = F,
		plotRegion = F,
		statSelectIn = F,
		statsCtrlsIn = F,
		statsRegion = F
		)
	values <- reactiveValues(
		newRatioFactor = NULL,
		oldRatioFactor = NULL,
		availableRatios = NULL,
		i = NULL
	)

########################################################  COVERAGE TAB  ################################################################

#~ ~~~~~~~~~~~~~~~~~~~~~~  BED: Load bed  ~~~~~~~~~~~~~~~~~~~~~~~~

	observeEvent(input$bedFile, {
		if (!is.null(input$bedFile)) {

				##  direct action
			names$bed <- gsub("\\.bed$","",input$bedFile$name)
			DFs$bed <- loadBed(input$bedFile$datapath,name=names$bed)
			toggles$bedCtrlsIn <-F
			if (toggles$coverage == T) { toggles$coverage <- F }
			if (toggles$isExampleCoverage == T) { toggles$isExampleCoverage <- F }

				##  changes to the side panel
			removeUI(selector='#analyseButton')
			removeUI(selector='#resetBed')
			insertUI(selector='#analyseOrReset',where="afterBegin",ui=actionButton("analyseButton","Analyse"))

				##  changes to the main panel
			removeUI(selector='#coverageDescription')
			plots$coveragePlot <- NULL
			insertUI(selector='#coverageMain',where="afterBegin",ui=tableOutput("bedContent"))
		}
	})
		##  outputs
	output$bedContent <- renderTable(DFs$bed[1:7,])

#~ ~~~~~~~~~~~~~~~~~~~~  BED: Load example  ~~~~~~~~~~~~~~~~~~~~~~~

	observeEvent(input$exampleCoverage, {
		removeUI(selector='#coverageDescription')
		toggles$isExampleCoverage <- T
		if (toggles$bedCtrlsIn == F) toggles$bedCtrlsIn <- T
	})

#~ ~~~~~~~~~~~~~~~~~~~~  BED: Analyse bed  ~~~~~~~~~~~~~~~~~~~~~~~~

	observeEvent(input$analyseButton, {
		if (toggles$bedCtrlsIn == F) toggles$bedCtrlsIn <- T
	})

#~ ~~~~~~~~~~~~~~~~~~  BED: Interface toggle  ~~~~~~~~~~~~~~~~~~~~~

	observeEvent(toggles$bedCtrlsIn, {
		if (toggles$bedCtrlsIn == T) {
			if (toggles$coverage == T) { toggles$coverage <- F }
			if (toggles$isExampleCoverage == T) {
				DFs$bed <- W303_S
				names$bed <- "W303_S"
			}
				##  changes to the side panel
			removeUI(selector='#analyseButton')
			removeUI(selector='#exampleCoverageDiv')
			insertUI(selector='#coverageSide',where="beforeEnd",ui=HTML(
				"<div id='bedSideCtrls'>

					<div id='rmChr'>
						<hr>
						<div class='myTooltip'><label>Remove chromosome</label><span class='myTooltiptext'>
							This allows you to remove all data from individual chromosomes
							</span></div>
						<div class='description'>Mitochondrial DNA is usually excluded from analysis.</div>
						<div style='padding-left:10px;'>
							<div id='rmChrInputDiv' class='inline' style='width:30%;margin-bottom:-5px;min-width:75px;'></div>
							<div id='rmChrButtonDiv' class='inline'></div>
						</div>
					</div>
					<div id='rmMax'>
						<hr>
						<div class='myTooltip'><label>Remove max value</label><span class='myTooltiptext'>
							Use this to remove one or more top outliers
							</span></div>
						<div class='description'>Remove individual outlier(s) with highest score.</div>
						<div style='padding-left:10px;'>
							<div id='rmMaxTimesInput' class='inline' style='width:20%;min-width:50px;'></div>
							<div id='rmMaxButton' class='inline'></div>
						</div>
					</div>
					<div id='rmOutliersDiv'>
						<hr>
						<div class='myTooltip'><label>Remove outliers (IQR)</label><span class='myTooltiptext'>
							Only use this if data is noisy for all the chromosomes
							</span></div>
						<div class='description'>
							Outliers (highlighted in grey) are either 3×IQR (interquartile range) or more above the third quartile
							or 3×IQR or more below the first quartile.<br>
							<b>Do not use</b> if only a few chromosomes appear noisy.<br>
							<b>Do not use</b> more than once.
						</div>
						<div id='rmOutliers' style='padding-left:10px;'></div>
					</div>
					<div id='rmOutliersManDiv'>
						<hr>
							<div class='myTooltip'><label>Remove outliers (median)</label><span class='myTooltiptext'>
							Typically, values less than 0.25*<b>median</b> are removed
							</span></div>
						<div class='description'>
							Specify minimum and/or maximum proportion to dataset median (pink line across) each bin should have.
							Bins that contain fewer reads per bin than Min*Median or more that Max*Median will be discarded.
						</div>
						<div id='rmOutliersManMin' class='inline' style='padding-left:10px;width:20%;'></div>
						<div id='rmOutliersManMax' class='inline' style='padding-left:5px;width:20%;'></div>
						<div id='rmOutliersManBtn' class='inline' style='padding-left:20px;width:40%;'></div>
					</div>



					<div>
						<div id='scatterPlotButtonDiv' style='text-align: center;'><hr></div>
					</div>
				</div>"
			))

			insertUI(selector='#rmChrInputDiv',where="beforeEnd",ui=uiOutput("chrRmOutput"))
			insertUI(selector='#rmChrButtonDiv',where="beforeEnd",ui=actionButton("rmChrButton","Remove chromosome"))
			insertUI(selector='#rmMaxTimesInput',where="beforeEnd",ui=textInput("rmMaxTimes",NULL,1,width="90%"))
			insertUI(selector='#rmMaxButton',where="beforeEnd",ui=actionButton("rmMaxButton","Remove max values"))
			insertUI(selector='#rmOutliers',where="beforeEnd",ui=actionButton("rmOutButton","Remove outliers"))
			insertUI(selector='#rmOutliersManMin',where="afterBegin",ui=textInput(
				"rmOutliersManMinInput",NULL,NULL,width="90%",placeholder="Min"
			))
			insertUI(selector='#rmOutliersManMax',where="afterBegin",ui=textInput(
				"rmOutliersManMaxInput",NULL,NULL,width="90%",placeholder="Max"
			))
			insertUI(selector='#rmOutliersManBtn',where="afterBegin",ui=actionButton("rmOutliersManButton","Remove outliers"))
			insertUI(
				selector='#scatterPlotButtonDiv',
				where="beforeEnd",
				ui=actionButton("scatterPlotButton","Switch to scatter plot view")
			)

				##  changes to the main panel
			removeUI(selector='#bedContent')
			insertUI(selector='#coverageMain',where="afterBegin",ui=HTML("
				<div id='bedMainCtrls'>
					<div id='coveragePlotDiv'></div>
					<div id='outlierInfoDiv' style='text-align:center;padding-left:5%;max-width:1000px;'>
						<div id='outlierDiv'><label>Hover over datapoints to display their properties</label></div>
					</div>
				</div>
			"))

			insertUI(
				selector='#coveragePlotDiv',
				where="afterBegin",
				ui=plotOutput(
					"coveragePlot",
					hover=hoverOpts(id="outliers"),
					width="100%",
					height=500
				)
			)
			insertUI(selector='#outlierDiv',where="beforeEnd",ui=verbatimTextOutput("outlierInfo"))

			insertUI(selector='#bedMainCtrls',where="beforeEnd",ui=HTML(
				"<div id='saveBedDiv'>
					<div id='saveBedFurther' class='inline''>
						<b style='display:block;padding: 20px;'>Save the data for further analysis:</b>
						<div id='bedSaveSelect' class='inline' style='width:20%;min-width:100px;'></div>
						<div id='bedSaveName' class='inline' style='width:30%;min-width:100px;'></div>
						<div id='bedSaveButton' class='inline'></div>
					</div>
					<div class='inline'>
						<b style='display:block;padding: 20px;'>Save data locally:</b>
						<div id='downloadCoveragePlot' class='inline'></div>
					</div>
				</div>"))
			insertUI(
				selector='#bedSaveSelect',
				where="afterBegin",
				ui=selectInput("bedType",NULL,c("Sample type"="","Replicating"="rep","Non-replicating"="nonRep"))
			)

			insertUI(selector='#bedSaveName',where="afterBegin",ui=textInput("bedName",NULL,names$bed))

			insertUI(selector='#bedSaveButton',where="afterBegin",ui=actionButton("saveBed","Save data"))

			insertUI(selector='#downloadCoveragePlot',where="afterBegin",ui=downloadButton("downloadCoveragePlot","Download plot"))
			insertUI(selector='#downloadCoveragePlot',where="beforeEnd",ui=downloadButton("downloadCoverageCSV","Download data"))

			names$plot <- names$bed
				##  direct action
			plots$coveragePlot <- plotBed(DFs$bed,plotting=F)
		} else {
			removeUI(selector='#bedSideCtrls')
			removeUI(selector='#bedMainCtrls')
		}
	})
		##  outputs
	output$coveragePlot <- renderPlot(plots$coveragePlot)
	output$outlierInfo <- renderPrint({
		nearPoints(DFs$bed,input$outliers,threshold=10,maxpoints=1)
	})

#~ ~~~~~~~~~~~~~~~~  BED: Cleanse bed  ~~~~~~~~~~~~~~~~~~

	observeEvent(input$rmChrButton, {
		if (input$chrRm!="") {
			DFs$bed <- rmChr(DFs$bed,input$chrRm)
			if (toggles$bedEdited == F) {
				insertUI(selector='#analyseOrReset',where="afterBegin",ui=actionButton("resetBed","Reset"))
				toggles$bedEdited <- T
			}
			if (toggles$coverage == F) {
				plots$coveragePlot <- plotBed(DFs$bed)
			} else {
				plots$coveragePlot <- plotCoverage(DFs$bed)
			}
		}
	})
		##  outputs
	output$chrRmOutput <- renderUI({
		chrList <- levels(DFs$bed$chrom)
		selectInput("chrRm",NULL,c("Select"="",chrList),width="90%")
	})


	observeEvent(input$rmMaxButton, {
		if (input$rmMaxTimes!="") {
			DFs$bed <- rmOutliers(DFs$bed,"max",n=input$rmMaxTimes)
			if (toggles$bedEdited == F) {
				insertUI(selector='#analyseOrReset',where="afterBegin",ui=actionButton("resetBed","Reset"))
				toggles$bedEdited <- T
			}
			if (toggles$coverage == F) {
				plots$coveragePlot <- plotBed(DFs$bed)
			} else {
				plots$coveragePlot <- plotCoverage(DFs$bed)
			}
		}
	})

	observeEvent(input$rmOutButton, {
		DFs$bed <- rmOutliers(DFs$bed,"IQR",range=3)
		if (toggles$bedEdited == F) {
			insertUI(selector='#analyseOrReset',where="afterBegin",ui=actionButton("resetBed","Reset"))
			toggles$bedEdited <- T
		}
		if (toggles$coverage == F) {
			plots$coveragePlot <- plotBed(DFs$bed)
		} else {
			plots$coveragePlot <- plotCoverage(DFs$bed)
		}
	})


	observeEvent(input$rmOutliersManButton, {
		if (input$rmOutliersManMinInput!="" & input$rmOutliersManMaxInput=="") {
			DFs$bed <- rmOutliers(DFs$bed,"median",loLim=input$rmOutliersManMinInput)
		}
		if (input$rmOutliersManMinInput=="" & input$rmOutliersManMaxInput!="") {
			DFs$bed <- rmOutliers(DFs$bed,"median",hiLim=input$rmOutliersManMaxInput)
		}
		if (input$rmOutliersManMinInput!="" & input$rmOutliersManMaxInput!="") {
			DFs$bed <- rmOutliers(DFs$bed,"median",loLim=input$rmOutliersManMinInput,hiLim=input$rmOutliersManMaxInput)
		}
		if (toggles$bedEdited == F) {
			insertUI(selector='#analyseOrReset',where="afterBegin",ui=actionButton("resetBed","Reset"))
			toggles$bedEdited <- T
		}
		if (toggles$coverage == F) {
			plots$coveragePlot <- plotBed(DFs$bed)
		} else {
			plots$coveragePlot <- plotCoverage(DFs$bed)
		}
	})




#~ ~~~~~~~~~~~~~~~~  BED: Box/Scatter plot switch  ~~~~~~~~~~~~~~~~~~
	observeEvent(input$scatterPlotButton, {
			##  changes to the side panel
		removeUI(selector='#scatterPlotButton')
		removeUI(selector='#rmMax')
		removeUI(selector='#rmOutliersDiv')
		insertUI(
			selector="#scatterPlotButtonDiv",
			where="beforeBegin",
			ui=HTML(
				"<div id='coverageRegion'>
					<hr>
					<div class='myTooltip'><label>Select region to plot</label><span class='myTooltiptext'>
					Focus on a chromosome or a smaller reagion
					</span></div>
					<div  style='padding-left:10px;'>
						<div id='coverageChrDiv' class='inline' style='width:20%;'></div>
						<div id='coverageChrStartDiv' class='inline' style='width:25%;'></div>
						<div id='coverageChrEndDiv' class='inline' style='width:25%;'></div>
						<div id='coverageRegionButtonDiv' class='inline'></div>
					</div>
					<div style='padding-left:10px;'>
						<div id='resetRegion' class='inline'></div>
					</div>
				</div>"
			)
		)
		insertUI(selector="#coverageChrDiv",where="afterBegin",ui=uiOutput('coverageChrOut',width="90%"))
		insertUI(
			selector="#coverageChrStartDiv",
			where="afterBegin",
			ui=textInput('coverageChrStart',NULL,NULL,width="95%",placeholder='Start')
		)
		insertUI(
			selector="#coverageChrEndDiv",
			where="afterBegin",
			ui=textInput('coverageChrEnd',NULL,NULL,width="95%",placeholder='End')
		)
		insertUI(selector="#coverageRegionButtonDiv",where="afterBegin",ui=actionButton('coverageRegionButton',"Plot region"))
		insertUI(
			selector="#scatterPlotButtonDiv",
			where="beforeEnd",
			ui=actionButton('boxPlotButton',"Switch to box plot view")
		)

			##  changes to the main panel
		removeUI(selector='#coveragePlot')
		removeUI(selector='#outlierDiv')
		insertUI(
			selector='#bedMainCtrls',
			where="afterBegin",
			ui=plotOutput(
				"coveragePlot",
				width=990,
				height=1200
			)
		)
			##  direct action
		plots$coveragePlot <- plotCoverage(DFs$bed,plotting=F)
		toggles$coverage <- T
	})

	observeEvent(input$boxPlotButton, {
			##  changes to the side panel
		removeUI(selector='#coverageRegion')
		removeUI(selector='#boxPlotButton')
		insertUI(selector="#rmChr", where="afterEnd", ui=HTML("

			<div id='rmMax'>
				<hr>
				<div class='myTooltip'><label>Remove max value</label><span class='myTooltiptext'>
					Use this to remove one or more top outliers
					</span></div>
				<div class='description'>Remove individual outlier(s) with highest score.</div>
				<div style='padding-left:10px;'>
					<div id='rmMaxTimesInput' class='inline' style='width:20%;min-width:50px;'></div>
					<div id='rmMaxButton' class='inline'></div>
				</div>
			</div>
			<div id='rmOutliersDiv'>
				<hr>
				<div class='myTooltip'><label>Remove outliers (IQR)</label><span class='myTooltiptext'>
					Only use this if data is noisy for all the chromosomes
					</span></div>
				<div class='description'>
					Outliers (highlighted in grey) are either 3×IQR (interquartile range) or more above the third quartile
					or 3×IQR or more below the first quartile.<br>
					<b>Do not use</b> if only a few chromosomes appear noisy.<br>
					<b>Do not use</b> more than once.
				</div>
				<div id='rmOutliers' style='padding-left:10px;'></div>
			</div>
		"))
		insertUI(selector='#rmMaxTimesInput',where="beforeEnd",ui=textInput("rmMaxTimes",NULL,1,width="90%"))
		insertUI(selector='#rmMaxButton',where="beforeEnd",ui=actionButton("rmMaxButton","Remove max values"))

		insertUI(selector='#rmOutliers',where="beforeEnd",ui=actionButton("rmOutButton","Remove outliers"))

		insertUI(
				selector='#scatterPlotButtonDiv',
				where="beforeEnd",
				ui=actionButton("scatterPlotButton","Switch to scatter plot view")
		)
			##  changes to the main panel
		removeUI(selector='#coveragePlot')
		insertUI(
				selector='#coveragePlotDiv',
				where="afterBegin",
				ui=plotOutput(
					"coveragePlot",
					hover=hoverOpts(id="outliers"),
					width="100%",
					height=500
				)
		)
		insertUI(selector='#outlierInfoDiv',where="afterBegin",ui=HTML("
			<div id='outlierDiv'><label>Hoover over datapoints to display their properties</label></div>
		"))
		insertUI(selector='#outlierDiv',where="beforeEnd",ui=verbatimTextOutput("outlierInfo"))
			##  direct action
		plots$coveragePlot <- plotBed(DFs$bed,plotting=F)
		toggles$coverage <- F
		toggles$coverageRegion <- F
	})

#~ ~~~~~~~~~~~~~~~~  Region plotting  ~~~~~~~~~~~~~~~~~~
	observeEvent(input$coverageChrIn, {
		if (input$coverageChrIn!="") {
			chrom <- input$coverageChrIn
			chromEnd <- max(DFs$bed$chromEnd[DFs$bed$chrom == chrom])
			updateTextInput(session,'coverageChrStart',value=0)
			updateTextInput(session,'coverageChrEnd',value=chromEnd)
		}
	})

	observeEvent(input$coverageRegionButton, {
		req(input$coverageChrIn,input$coverageChrStart,input$coverageChrEnd)
		region <- paste0(input$coverageChrIn,":",input$coverageChrStart,"-",input$coverageChrEnd)
		if (toggles$coverageRegion==F) {
			removeUI(selector='#coveragePlot')
			insertUI(
				selector='#bedMainCtrls',
				where="afterBegin",
				ui=plotOutput(
					"coveragePlot",
					width="70%",
					height=300
				)
			)
			insertUI(selector="#resetRegion",where="afterBegin",ui=actionButton('resetRegionButton',"Reset view"))
		}
		plots$coveragePlot <- plotCoverage(DFs$bed,region=region,plotting=F)
		toggles$coverageRegion <- T
	})

	observeEvent(input$resetRegionButton, {
		removeUI(selector='#resetRegionButton')
		removeUI(selector='#coveragePlot')
		insertUI(
			selector='#bedMainCtrls',
			where="afterBegin",
			ui=plotOutput(
				"coveragePlot",
				width=990,
				height=1200
			)
		)
			##  direct action
		plots$coveragePlot <- plotCoverage(DFs$bed,plotting=F)
		toggles$coverageRegion <- F
	})
		##  outputs
	output$coverageChrOut <- renderUI({
		chrList <- levels(DFs$bed$chrom)
		selectInput("coverageChrIn",NULL,c("Chr"="",chrList))
	})

#~ ~~~~~~~~~~~~~~~~  Reset bed  ~~~~~~~~~~~~~~~~~~

	observeEvent(input$resetBed, {
			##  changes to the side panel
		removeUI(selector='#resetBed')
			##  update reactive values
		toggles$bedEdited <- F
			##  direct action
		if (toggles$isExampleCoverage == T) {
			names$bed <- "W303_S"
			DFs$bed <- W303_S
		} else {
			DFs$bed <- loadBed(input$bedFile$datapath,name=names$bed)
		}
		if (toggles$coverage == F ) {
			plots$coveragePlot <- plotBed(DFs$bed,plotting=F)
		} else {
			plots$coveragePlot <- plotCoverage(DFs$bed,plotting=F)
		}
	})

#~ ~~~~~~~~~~~~~~~~  Save bed  ~~~~~~~~~~~~~~~~~~

	observeEvent(input$saveBed, {
		req(input$bedType,input$bedName)
		names$bed <- input$bedName
		sampleName <- input$bedName
		DFs$bed$name <- rep(sampleName,dim(DFs$bed)[1])
		if (input$bedType == 'rep') {
			DFs$rep[[sampleName]] <- DFs$bed
				## add inteface to the ratio tab if we have one of each, rep and non-rep
			if (length(names(DFs$rep))==1 & length(names(DFs$nonRep)) >= 1) {
				insertUI(selector='#ratioSide',where="afterBegin",ui=HTML("<div id='ratioSelectorDiv' style='padding-bottom:20px;'></div>"))
				insertUI(selector="#ratioSelectorDiv",where="afterBegin",ui=uiOutput('nonRepSamples'))
				insertUI(selector="#ratioSelectorDiv",where="afterBegin",ui=uiOutput('repSamples'))
			}
		}
		if (input$bedType == 'nonRep') {
			DFs$nonRep[[sampleName]] <- DFs$bed
				## add inteface to the ratio tab if we have one of each, rep and non-rep
			if (length(names(DFs$nonRep))==1 & length(names(DFs$rep)) >= 1) {
				insertUI(selector='#ratioSide',where="afterBegin",ui=HTML("<div id='ratioSelectorDiv' style='padding-bottom:20px;'></div>"))
				insertUI(selector="#ratioSelectorDiv",where="afterBegin",ui=uiOutput('nonRepSamples'))
				insertUI(selector="#ratioSelectorDiv",where="afterBegin",ui=uiOutput('repSamples'))
			}
		}
		removeUI(selector='#bedSideCtrls')
		removeUI(selector='#resetBed')
		removeUI(selector='#saveBedFurther')
		##  update reactive values
		toggles$bedEdited <- F
		toggles$coverage <- F
		toggles$coverageRegion <- F
	})
		##  outputs
	output$downloadCoveragePlot <- downloadHandler(
		filename = function() {
			if (toggles$coverageRegion==F) {
				paste0(names$bed, '_RawReads.pdf')
			} else {
				region <- paste0(input$coverageChrIn,"_",input$coverageChrStart,"-",input$coverageChrEnd)
				paste0(names$bed,'_',region,'_RawReads.pdf')
			}
		},
		content = function(file) {
			if (toggles$coverage == F) {
				ggsave(file, plot=plots$coveragePlot, device=cairo_pdf, width = 40, height = 20, units = "cm")
			} else {
				if (toggles$coverageRegion==F) {
					ggsave(file, plot=plots$coveragePlot, device=cairo_pdf, width = 25, height = 35, units = "cm")
				} else {
					ggsave(file, plot=plots$coveragePlot, device=cairo_pdf, width = 30, height = 12, units = "cm")
				}
			}
		}
	)

	output$downloadCoverageCSV <- downloadHandler(
	    filename = function() {
			if (toggles$coverageRegion==F) {
				paste0(names$bed, '_RawReads.tsv')
			} else {
				region <- paste0(input$coverageChrIn,"_",input$coverageChrStart,"-",input$coverageChrEnd)
				paste0(names$bed,'_',region,'_RawReads.tsv')
			}
		},
	    content = function(file) {
			fileContent <- DFs$bed
			fileContent$name <- rep(input$bedName,dim(DFs$bed)[1])
			if (toggles$coverageRegion==F) {
				write.table(fileContent,file=file,sep="\t",col.names=T,row.names=F,quote=F)
			} else {
				fileContent <- subset(
					fileContent,
					chrom==input$coverageChrIn & chromStart>=input$coverageChrStart & chromEnd<=input$coverageChrEnd
					)
				write.table(fileContent,file=file,sep="\t",col.names=T,row.names=F,quote=F)
			}
	    })







#######################################################  RATIO TAB  ###################################################################






	observeEvent(input$repBed, {
		if (input$nonRepBed!="" && input$repBed!="") {
			removeUI(selector='#makeRatioButton')
			insertUI(selector="#nonRepSamples",where="afterEnd",ui=actionButton('makeRatioButton',"Make ratio"))
		}
	})
	output$repSamples <- renderUI({
		repSampleList <- names(DFs$rep)
		selectInput("repBed","Select samples to calculate ratio:",c("Replicating"="",repSampleList),width="75%")
	})


	observeEvent(input$nonRepBed, {
		if (input$repBed!="" && input$nonRepBed!="") {
			removeUI(selector='#makeRatioButton')
			insertUI(selector="#nonRepSamples",where="afterEnd",ui=actionButton('makeRatioButton',"Make ratio"))
		}
	})
	output$nonRepSamples <- renderUI({
		nonRepSampleList <- names(DFs$nonRep)
		selectInput("nonRepBed",NULL,c("Non-replicating"="",nonRepSampleList),width="75%")
	})

	observeEvent(input$ratioFile, {
		if (!is.null(input$ratioFile)) {
			DFs$ratio <- read.table(input$ratioFile$datapath,sep=",",header=T,colClasses=
				c("factor","integer","integer","factor","factor","numeric","factor"))
			plots$ratioPlot <- plotRatio(DFs$ratio$ratio,plotting=F)
			if (toggles$isExampleRatio == T) {
				toggles$isExampleRatio <- F
				insertUI('#saveRatioFurther',"afterBegin",ui=HTML("<b style='display:block;padding: 20px;'>Save the ratio for plotting:</b>"))
				insertUI('#saveRatioBut',"afterBegin",ui=actionButton("saveRatioButton","Save ratio"))
			}
			if (toggles$ratioCtrlsIn == T) {
				updateTextInput(session,"loLimInput",value="")
				updateTextInput(session,"hiLimInput",value="")
				values$newRatioFactor <- as.numeric(DFs$ratio$ratioFactor[1])
				updateTextInput(session,"ratioFactor",value=as.character(DFs$ratio$ratioFactor[1]))
			} else {
				toggles$ratioCtrlsIn <- T
			}
		}
	})

	observeEvent(input$exampleRatio, {
		toggles$isExampleRatio <- T
		toggles$ratioCtrlsIn <- T
	})

	observeEvent(input$makeRatioButton, {
		req(input$nonRepBed,input$repBed)
		removeUI(selector='#makeRatioButton')
		DFs$ratio <- makeRatio(DFs$rep[[input$repBed]],DFs$nonRep[[input$nonRepBed]])
		plots$ratioPlot <- plotRatio(DFs$ratio$ratio,plotting=F)
		if (toggles$isExampleRatio == T) {
			toggles$isExampleRatio <- F
			insertUI('#saveRatioFurther',"afterBegin",ui=HTML("<b style='display:block;padding: 20px;'>Save the ratio for plotting:</b>"))
			insertUI('#saveRatioBut',"afterBegin",ui=actionButton("saveRatioButton","Save ratio"))

		}
		if (toggles$ratioCtrlsIn == T) {
			updateTextInput(session,"loLimInput",value="")
			updateTextInput(session,"hiLimInput",value="")
			values$newRatioFactor <- 1.000
			updateTextInput(session,"ratioFactor",value="1.000")
		} else {
			toggles$ratioCtrlsIn <- T
		}

	})

	observeEvent(DFs$ratio, {
		names$currentRatio <- c(as.character(DFs$ratio$name.rep[1]),as.character(DFs$ratio$name.nonRep[1]))
	},ignoreInit=T)

	observeEvent(toggles$ratioCtrlsIn, {
		if (toggles$ratioCtrlsIn == T) {
			removeUI('#ratioDescription')
			if (toggles$isExampleRatio == T) {
				DFs$ratio <- W303
			}
			ratioFactor <- as.character(DFs$ratio$ratioFactor[1])
				##  changes to the side panel
			removeUI(selector='#makeRatioButton')
			removeUI(selector='#exampleRatioDiv')
			insertUI(
				selector='#ratioSide',
				where="afterEnd",
				ui=HTML("
					<div id='ratioSideCtrls'>
						<div id='ratioTrimDiv'>
							<hr>
							<div class='myTooltip'><label>Trim the ratio</label><span class='myTooltiptext'>
							Use before automatic normalisation. 0.5-1.5 is a safe starting point
							</span></div>
							<div class='description'>
								Some genomic regions may exhibit high variability in sequencing depth.
								Exclude them by trimming the ratio values.
							</div>
							<div style='padding-left:10px;'>
								<div id='loLimDiv' class='inline' style='width:30%;'></div>
								<div id='hiLimDiv' class='inline' style='width:30%;'></div>
								<div id='trimButtonDiv' class='inline' style='padding-left:10px;'></div>
							</div>
						</div>
						<div id='autoNormDiv'>
							<hr>
							<div class='myTooltip'><label>Automatic normalisation</label><span class='myTooltiptext'>
							Use with trimmed full range S phase samples
							</span></div>
							<div class='description'>
								This fits the data on a scale from 1 to a value entered below (2 for sorted samples, between 1 and 3 for MFA-seq), by minimising the sum of the outliers.
							</div>
							<div style='padding-left:10px;'>
								<div id='upperLimitField' class='inline' style='width:30%;'></div>
								<div id='autoNormBtn' class='inline'></div>
							</div>
						</div>
						<div id='maNormDiv'>
							<hr>
							<div class='myTooltip'><label>Manual normalisation</label><span class='myTooltiptext'>
							Use with asynchronous or early S phase cell samples
							</span></div>
							<div class='description'>
								If the automatic normalisation was used, it will show the calculated value.
							</div>
							<div style='padding-left:10px;'>
								<div id='maNormField' class='inline' style='width:30%;'></div>
								<div id='maNormButton' class='inline'></div>
							</div>
						</div>
					</div>
				")
			)
			insertUI(selector="#loLimDiv",where="afterBegin",ui=textInput('loLimInput',NULL,placeholder="Low limit"))
			insertUI(selector="#hiLimDiv",where="afterBegin",ui=textInput('hiLimInput',NULL,placeholder="High limit"))
			insertUI(selector="#trimButtonDiv",where="afterBegin",ui=actionButton('trimRatioButton',"Trim"))

			insertUI(selector="#upperLimitField",where="afterBegin",ui=textInput('upperLimit',NULL,"2.0",width="90%"))

			insertUI(selector="#autoNormBtn",where="afterBegin",ui=actionButton('normaliseButton',"Auto normalise"))
			insertUI(selector="#maNormField",where="afterBegin",ui=textInput('ratioFactor',NULL,ratioFactor,width="90%"))
			insertUI(selector="#maNormButton",where="afterBegin",ui=actionButton('maNormButton',"Update"))
			values$newRatioFactor <- as.numeric(ratioFactor)

				##  changes to the main panel
			insertUI(
				selector='#ratioMain',
				where="afterBegin",
				ui=HTML(paste0("
					<div id='ratioMainCtrls'>
						<div id='ratioPlotDiv'></div>
						<div id='saveRatioDiv' style='padding: 1cm 0.5cm;'>
					<div id='saveRatioFurther' class='inline' style='padding-right:2cm;'>",
						if (toggles$isExampleRatio == F) {
							paste0("<b style='display:block;padding: 20px;'>Save the ratio for plotting:</b>") },
								"<div id='saveRatioBut'></div>
							</div>
							<div class='inline'><b style='display:block;padding: 20px;'>Save data locally:</b>
								<div id='downloadRatio' class='inline'></div>
							</div>
						</div>
					</div>
				"))
			)

			insertUI(selector='#ratioPlotDiv',where="afterBegin",ui=plotOutput("plotHist",width="80%",height=500))
			plots$ratioPlot <- plotRatio(DFs$ratio$ratio,plotting=F)
			if (toggles$isExampleRatio != T) {
				insertUI(selector='#saveRatioBut',where="afterBegin",ui=actionButton("saveRatioButton","Save ratio"))
			}
			insertUI(selector='#downloadRatio',where="afterBegin",ui=downloadButton("downloadRatioPlot","Download plot"))
			insertUI(selector='#downloadRatio',where="beforeEnd",ui=downloadButton("downloadRatioCSV","Download data"))
		} else {
			removeUI(selector='#ratioSideCtrls')
			removeUI(selector='#ratioMainCtrls')
		}
	})

	output$plotHist <- renderPlot(plots$ratioPlot)

	observeEvent(input$trimRatioButton, {
		req(input$loLimInput,input$hiLimInput)
		loLim <- as.numeric(input$loLimInput)
		hiLim <- as.numeric(input$hiLimInput)
		DFs$ratio <- trimRatio(DFs$ratio,loLim,hiLim)
		if ("tmpRatio" %in% colnames(DFs$ratio)) {
			plots$ratioPlot <- plotRatio(DFs$ratio$tmpRatio,plotting=F)
		} else {
			plots$ratioPlot <- plotRatio(DFs$ratio$ratio,plotting=F)
		}
	})

	observeEvent(input$normaliseButton, {
		DFs$ratio <- normaliseRatio(DFs$ratio,upperLimit=as.numeric(input$upperLimit),replace=F)
		plots$ratioPlot <- plotRatio(DFs$ratio$tmpRatio,plotting=F)
		values$oldRatioFactor <- values$newRatioFactor
		values$newRatioFactor <- as.numeric(attributes(DFs$ratio)$comment)
		updateTextInput(session,"ratioFactor",value=values$newRatioFactor)
		updateTextInput(session,"loLimInput",value=round(values$newRatioFactor*as.numeric(input$loLimInput)/values$oldRatioFactor,2))
		updateTextInput(session,"hiLimInput",value=round(values$newRatioFactor*as.numeric(input$hiLimInput)/values$oldRatioFactor,2))

	})

	observeEvent(input$maNormButton, {
		values$oldRatioFactor <- values$newRatioFactor
		values$newRatioFactor <- as.numeric(input$ratioFactor)
		DFs$ratio <- normaliseRatio(DFs$ratio,rFactor=input$ratioFactor,replace=F)
		plots$ratioPlot <- plotRatio(DFs$ratio$tmpRatio,plotting=F)
		updateTextInput(session,"loLimInput",value=round(values$newRatioFactor*as.numeric(input$loLimInput)/values$oldRatioFactor,2))
		updateTextInput(session,"hiLimInput",value=round(values$newRatioFactor*as.numeric(input$hiLimInput)/values$oldRatioFactor,2))
	})




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  changes to the PLOT tab  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


	observeEvent(input$examplePlot, {
		toggles$isExamplePlot <- T
		if (toggles$plotCtrlsIn == F) toggles$plotCtrlsIn <- T
		DFs$ratios <- syncSeq[["data"]]
		DFs$guide <- syncSeq[["guide"]]
	})

		##  Inserting plotting interface either after first ratio in or when using example plot
	observeEvent(toggles$plotCtrlsIn, {
		if (toggles$plotCtrlsIn == T) {
			removeUI(selector='#examplePlotDiv')
			insertUI(selector='#samples',where="afterBegin",ui=HTML("
				<div class='myTooltip' style='padding-bottom:15px;'><ctrlH>Samples</ctrlH>
					<span class='myTooltiptext'>This area controls appearance of individual samples</span>
				</div>
				<div id='sampleArea' style='padding-left:10px;'></div>
			"))
			insertUI(
				selector='#sampleArea',
				where="afterEnd",
				ui=HTML("
					<hr>
					<div id='plotElements'>
						<div class='myTooltip'><ctrlH>Additional features</ctrlH>
							<span class='myTooltiptext'>Add features to the plot using your data in bed format</span>
						</div>
						<div style='padding:15px 10px;'>
							<div id='selectElement' class='inline' style='width:35%;'></div>
							<div id='uploadElement' class='inline' style='width:60%;padding-left:15px;'></div>
							<div id='plotFeaturesDiv'></div>
						</div>
					</div>
					<hr>
					<div style='padding-bottom:15px;'>
						<div class='myTooltip'><ctrlH>Smoothing controls</ctrlH>
							<span class='myTooltiptext'>Noisy data may benefit from adding a smoother</span>
						</div>
						<div id='ratioSmoothControls' style='padding:15px 10px 15px 10px;'>
							<div id='groupDiv' class='inline' style='width:35%;'>
								<div class='inline myTooltip' style='padding-top:7px;'><label>Group size:</label>
									<span class='myTooltiptext'>Minimum number of bins in a group</span>
								</div>
								<div id='smoothGroup' class='inline' style='width:40%;'></div>
							</div>
							<div id='splitDiv' class='inline' style='width:30%;padding-left:10px;'>
								<div class='inline myTooltip' style='padding-top:7px;'><label>Split:</label>
									<span class='myTooltiptext'>Number of missing bins along a chromosome to initiate a new group</span>
								</div>
								<div id='smoothSplit' class='inline' style='width:40%;'></div>
							</div>
							<div id='smoothBtn' class='inline' style='width:25%;'></div>
						</div>
					</div>
					<hr>
					<div>
						<div class='myTooltip'><ctrlH>Plotting controls</ctrlH>
							<span class='myTooltiptext'>Choose the appearance of the plot</span>
						</div>
						<div id='plottingCtrls' style='padding:10px;'>
							<div id='plotRegionDiv' style='padding-bottom:10px;'>
								<div id='plotRegion'>
									<b>Select region to plot</b><br>
									<div>
										<div id='plotChrDiv' class='inline' style='width:20%;'></div>
										<div id='plotChrStartDiv' class='inline' style='width:23%;'></div>
										<div id='plotChrEndDiv' class='inline' style='width:23%;'></div>
										<div id='resetPlotRegionDiv' class='inline' style='width:28%;'></div>
									</div>
								</div>
							</div>
							<div id='geomInputDiv' class='inline' style='width:45%;'>
								<div style='padding-top:5px;'><b>Plot type:</b></div>
								<div id='geomDiv' style='padding-right:10px;width:80%;'></div>
							</div>
							<div id='plotLimits' class='inline' style='width:45%;'>
								<div style='padding-top:5px;'><b><i>y</i> axis limits:</b></div>
								<div>
									<div id='plotLim1' class='inline' style='width:45%;'></div>
									<div id='plotLim2' class='inline' style='width:45%;'></div>
								</div>
							</div>
							<div id='plotGenomeBtnDiv' style='padding-top:20px;text-align:center;'></div>
						</div>
					</div>
				")
			)

			insertUI(selector="#selectElement",where="afterBegin",ui=selectInput('plotFeatures',NULL,
				c("Choose type:"="","vLine"="lines","Circle"="circles","Rectangle"="rectangles","Pointer"="pointers")))
			insertUI(selector="#uploadElement",where="afterBegin",ui=fileInput(
				'plotFeaturesFile',NULL,multiple=F,accept=".bed",buttonLabel = "Browse...", placeholder="No file selected"
			))
			insertUI(selector="#smoothGroup",where="afterBegin",ui=textInput('group',NULL,4))
			insertUI(selector="#smoothSplit",where="afterBegin",ui=textInput('split',NULL,5))
			insertUI(selector="#smoothBtn",where="afterBegin",ui=actionButton('smoothButton',"Apply smoothing"))

			insertUI(selector="#plotChrDiv",where="afterBegin",ui=uiOutput('plotChrOut',width="90%"))
			insertUI(
				selector="#plotChrStartDiv",
				where="afterBegin",
				ui=textInput('plotChrStart',NULL,NULL,width="95%",placeholder='Start')
			)
			insertUI(
				selector="#plotChrEndDiv",
				where="afterBegin",
				ui=textInput('plotChrEnd',NULL,NULL,width="95%",placeholder='End')
			)
			insertUI(selector='#resetPlotRegionDiv',where="afterBegin",ui=actionButton("resetPlotRegion","Reset region"))
			if (toggles$isExamplePlot) {
				geom_default <- 'geom_segment'
			} else {
				geom_default <- 'geom_point'
			}
			insertUI(selector="#geomDiv",where="afterBegin",ui=selectInput(
				'geomInput',
				NULL,
				c("scatter"='geom_point',"polygon"='geom_ribbon',"bar"='geom_segment'),
				selected=geom_default)
			)
			insertUI(
				selector="#plotLim1",
				where="afterBegin",
				ui=textInput('plotLimLow',NULL,1,width="95%",placeholder='From')
			)
			insertUI(
				selector="#plotLim2",
				where="afterBegin",
				ui=textInput('plotLimHi',NULL,2,width="95%",placeholder='To')
			)
			insertUI(selector="#plotGenomeBtnDiv",where="afterBegin",ui=actionButton('plotGenomeButton'," Plot "))
		}
	})




	observeEvent(input$saveRatioButton, {
		if (toggles$isExampleRatio == T) toggles$isExampleRatio <- F
		if (toggles$plotCtrlsIn == F) toggles$plotCtrlsIn <- T

	##  direct action
		if ("tmpRatio" %in% colnames(DFs$ratio)) {
			DFs$ratio$ratio <- round(DFs$ratio$tmpRatio,3)
			DFs$ratio$tmpRatio <- NULL
		}
		names$ratio[nrow(names$ratio)+1,] <- c(DFs$ratio$name.rep[1],DFs$ratio$name.nonRep[1],paste0(DFs$ratio$name.rep[1]," (",DFs$ratio$name.nonRep[1],")"))
		DFs$ratio$group <- paste("NA")
		DFs$ratio$splineSmooth <- paste("NA")
		DFs$ratio$ratioFactor <- as.character(round(values$newRatioFactor,3))
		DFs$ratio$ratioFactor <- factor(DFs$ratio$ratioFactor,levels=unique(DFs$ratio$ratioFactor))
		DFs$ratios <- rbind(DFs$ratios,DFs$ratio)
		values$availableRatios <- unique(paste0(as.character(DFs$ratios$name.rep)," (",as.character(DFs$ratios$name.nonRep),")"))
		DFs$ratios2add <- unique(DFs$ratio[,c("name.rep","name.nonRep")])

	##  changes to the side panel
		toggles$ratioCtrlsIn <- F
	##  changes to the main panel
	})



	observeEvent(values$i, {
		if (values$i ==2) {
			if (toggles$statSelectIn == F) toggles$statSelectIn <- T
		}
	})


	output$downloadRatioPlot <- downloadHandler(
	    filename = function() { paste0(names$currentRatio[1],'_(',names$currentRatio[2],')_RatioHist.pdf') },
	    content = function(file) {
		        ggsave(file, plot=plots$ratioPlot, device=cairo_pdf, width = 30, height = 15, units = "cm")
	    }
	)
	output$downloadRatioCSV <- downloadHandler(
	    filename = function() {	paste0(names$currentRatio[1],'_(',names$currentRatio[2],')_Ratio.csv') },
	    content = function(file) {
			ratioData <- DFs$ratio
			ratioData$ratioFactor <- as.character(round(values$newRatioFactor,3))
			ratioData$ratioFactor <- factor(ratioData$ratioFactor,levels=unique(ratioData$ratioFactor))
			if ("tmpRatio" %in% colnames(ratioData)) {
				ratioData$ratio <- round(ratioData$tmpRatio,3)
				ratioData$tmpRatio <- NULL
			}
			write.table(ratioData,file=file,sep=",",col.names=T,row.names=F,quote=F)
	    }
	)




###########################################################~~~~~~~~~~~~################################################################
###########################################################  PLOT TAB  ################################################################
###########################################################~~~~~~~~~~~~################################################################




	observeEvent(input$plotFile, {
		if (!is.null(input$plotFile)) {
			if (toggles$isExamplePlot == T) toggles$isExamplePlot <- F
			if (toggles$plotCtrlsIn == F) {	toggles$plotCtrlsIn <- T }
			plotFileContent <- read.table(input$plotFile$datapath,sep=",",header=T,colClasses=c("factor","integer","integer","factor","factor","numeric","factor","integer","numeric"))
			samples <- unique(plotFileContent[,c("name.rep","name.nonRep")])
			samples$name <- paste0(as.character(samples[,1])," (",as.character(samples[,2]),")")
			names$ratio <- rbind(names$ratio,samples)
			DFs$ratios2add <- unique(plotFileContent[,c("name.rep","name.nonRep")])
			DFs$ratios <- rbind(DFs$ratios,plotFileContent)
		}
	})

	observeEvent(input$examplePlot, {
		toggles$isExamplePlot <- T
		DFs$ratios <- syncSeq[["data"]]
		if (toggles$plotCtrlsIn == F) toggles$plotCtrlsIn <- T
		DFs$ratios2add <- unique(DFs$ratios[,c("name.rep","name.nonRep")])
		plots$genomePlot <- plotGenome(DFs$ratios,plotting=F,guide=syncSeq[["guide"]],geom="geom_segment")
			##  changes to main panel
		removeUI(selector='#plotGenomeOut')
		removeUI(selector='#downloadButtons')
		removeUI(selector='#plotDescription')
		insertUI(selector='#plotMain',where="afterBegin",ui=plotOutput('plotGenomeOut',width=990,height=1200))
		insertUI(selector='#plotGenomeOut',where="afterEnd",ui=tags$div(id = 'downloadButtons'))
		insertUI(selector='#downloadButtons',ui=downloadButton("downloadPlot","Download plot"))
		insertUI(selector='#downloadButtons',ui=downloadButton("downloadCSV","Download data"))
	})


	observeEvent(DFs$ratios2add, {
		if (nrow(DFs$ratios2add) != 0) {
			if (is.null(values$i)) values$i <- as.integer(0)
			for (i in 1:nrow(DFs$ratios2add)) {
					##  generate/add to the guide dataframe for plotting
				colors = c("#7EA4D6","#B7234B","#1F00F0","#F8D733","#00811F","#E680E6","#808080")
				if (!toggles$isExamplePlot) {
					DFs$guide[as.integer(values$i)+i,] <- c(
						as.integer(as.integer(values$i)+i),
						as.character(DFs$ratios2add$name.rep[i]),
						as.character(DFs$ratios2add$name.nonRep[i]),
						as.logical(TRUE),
						as.logical(FALSE),
						as.character(colors[as.integer(values$i)+i])
					)
				}
			##  changes to the PLOT tab
				insertUI(selector='#sampleArea',where="beforeEnd",ui=HTML(paste0("
					<div id='ratioDiv",as.integer(values$i)+i,"' style='color:#696969;padding-bottom:10px;'>
						<b
							style='display:block;padding-bottom:10px;color:#000000;'>
							<i>",as.character(DFs$ratios2add$name.rep[i])," (",as.character(DFs$ratios2add$name.nonRep[i]),")</i>
						</b>","
						<div id='orderDiv",as.integer(values$i)+i,"' class='inline' style='width:15%;'>
							<div id='tmpOrderDiv",as.integer(values$i)+i,"'></div>
						</div>
						<div class='inline'><b style='display:block;padding-left:10px;'>Data</b>
							<div class='inline'>
								<div id='rawDiv",as.integer(values$i)+i,"' class='inline' style='padding-left:10px;'></div>
								<div class='inline' style='display:inline-block;padding-top:12px;'>Raw</div>
							</div>
							<div class='inline'>
								<div id='smoothDiv",as.integer(values$i)+i,"' class='inline' style='padding-left:12px;'></div>
								<div id='smoothie",as.integer(values$i)+i,"' class='inline' style='display:inline-block;padding-top:12px;'>
									<div id='rmSmooth",as.integer(values$i)+i,"'",
										if (toggles$isExamplePlot!=T) {paste0(" style='text-decoration:line-through;'")},
										">Smooth
									</div>
								</div>
							</div>
						</div>
						<div class='inline'><b style='display:block;padding-left:10px;'>Colour</b>
							<div id='colorDiv",as.integer(values$i)+i,"' style='padding-top:5px;padding-left:10px;'> </div>
						</div>

					</div>")
					)
				)
				insertUI(selector=paste0('#tmpOrderDiv',as.integer(values$i)+i),where="afterBegin",ui=selectInput(
					paste0('orderInput',as.integer(values$i)+i),
					"Order",
					c(seq(1,as.integer(values$i)+i),NA),
					selected=as.integer(values$i)+i)
				)
				insertUI(selector=paste0('#rawDiv',as.integer(values$i)+i),where="afterBegin",ui=checkboxInput(
					paste0('rawInput',as.integer(values$i)+i),
					NULL,
					TRUE)
				)
				insertUI(selector=paste0('#smoothDiv',as.integer(values$i)+i),where="afterBegin",ui=checkboxInput(
					paste0('smoothInput',as.integer(values$i)+i),
					NULL,
					FALSE)
				)
				insertUI(selector=paste0('#colorDiv',as.integer(values$i)+i),where="afterBegin",ui=colourpicker::colourInput(
					paste0('colorInput',as.integer(values$i)+i),
					NULL,
					paste(colors[as.integer(values$i)+i])
					)
				)
			}
			if (as.integer(values$i)+nrow(DFs$ratios2add) > 1) {
				k <- (as.integer(as.integer(values$i)+nrow(DFs$ratios2add)))
				for (j in 1:(k-1)) {
					removeUI(paste0("#tmpOrderDiv",j))
					insertUI(paste0("#orderDiv",j),"afterBegin",HTML(paste0("<div id='tmpOrderDiv",j,"'></div>")))
					insertUI(paste0("#tmpOrderDiv",j),"afterBegin",
						selectInput(paste0('orderInput',j),"Order",c(seq(1,k),NA),selected=j))
				}
			}
			if (toggles$isExamplePlot != T) {
				if (as.integer(values$i) < 2 & as.integer(values$i) + i > 2) {
					if (toggles$statSelectIn == F) toggles$statSelectIn <- T
				}
				values$i <- as.integer(values$i) + i
			}
		}
	})

	observeEvent(input$plotFeaturesFile, {
		if (input$plotFeatures!="") {
			filePath <- input$plotFeaturesFile$datapath
			feature_str <- paste0("DFs$",input$plotFeatures,"<-loadBed(filePath)")
			eval(parse(text=feature_str))
			str <- paste0("featureName <- levels(DFs$",input$plotFeatures,"$name)[1]")
			eval(parse(text=str))
			features <- c("lines","circles","rectangles","pointers")
			colours <- c("#00FF00","#FFFFFF","#FF0000","#FF7F00")
			removeUI(paste0("#",input$plotFeatures))
			insertUI('#plotFeaturesDiv',"beforeEnd",ui=HTML(paste0("
				<div class='inline'>
					<div id='",input$plotFeatures,"' class='inline' style='max-width:75%;padding-left:10px;'></div>
					<div id='colour",input$plotFeatures,"' class='inline' style='max-width:25%;padding-left:10px;'></div>
				</div>"
			)))
			insertUI(paste0("#",input$plotFeatures),"afterBegin",checkboxInput(paste(input$plotFeatures),paste(featureName),TRUE))
			insertUI(paste0("#colour",input$plotFeatures),"afterBegin",colourpicker::colourInput(
				paste0(input$plotFeatures,"ColourInput"),NULL,colours[match(input$plotFeatures,features)],showColour = "background",palette = "limited"))

		}
	})

	observeEvent(input$plotChrIn, {
		if (input$plotChrIn!="") {
			chrom <- input$plotChrIn
			chromEnd <- max(DFs$ratios$chromEnd[DFs$ratios$chrom == chrom])
			updateTextInput(session,'plotChrStart',value=0)
			updateTextInput(session,'plotChrEnd',value=chromEnd)
			toggles$plotRegion <- T
		}
	})

	observeEvent(input$resetPlotRegion, {
		updateSelectInput(session,'plotChrIn',selected="")
		updateTextInput(session,'plotChrStart',value="")
		updateTextInput(session,'plotChrEnd',value="")
		toggles$plotRegion <- F
	})

	observeEvent(input$plotGenomeButton, {
		req(DFs$ratios,DFs$guide,input$plotLimLow,input$plotLimHi)
		removeUI(selector='#plotGenomeOut')
		removeUI(selector='#downloadButtons')
		removeUI(selector='#plotDescription')
	##  direct action: read inputs
		DFs$guide$order <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.integer(input[[paste0("orderInput",i)]])
		})) #order
		DFs$guide$raw <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.logical(input[[paste0("rawInput",i)]])
		})) #raw
		DFs$guide$smooth <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.logical(input[[paste0("smoothInput",i)]])
		})) #smooth
		DFs$guide$color <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.character(input[[paste0("colorInput",i)]])
		})) #color
		ylims <- c(as.numeric(input$plotLimLow),as.numeric(input$plotLimHi))

	##  direct action: create plot object
			##  check if region is selected
		plotString <- "plots$genomePlot <- plotGenome(DFs$ratios,plotting=F,guide=DFs$guide,geom=input$geomInput,ylims=ylims"
		if (input$plotChrIn!="") {
			plotString <- paste0(plotString,",region='",input$plotChrIn,":",input$plotChrStart,"-",input$plotChrEnd,"'")
		}
		if (!is.null(DFs$lines)) {
			if (input$lines) {
				plotString <- paste0(plotString,",lines=DFs$lines,colourLines=input$linesColourInput")
			}
		}
		if (!is.null(DFs$circles)) {
			if (input$circles) {
				plotString <- paste0(plotString,",circles=DFs$circles,colourCircles=input$circlesColourInput")
			}
		}
		if (!is.null(DFs$rectangles)) {
			if (input$rectangles) {
				plotString <- paste0(plotString,",rectangles=DFs$rectangles,colourRectangles=input$rectanglesColourInput")
			}
		}
		if (!is.null(DFs$pointers)) {
			if (input$pointers) {
				plotString <- paste0(plotString,",pointers=DFs$pointers,colourPointers=input$pointersColourInput")
			}
		}
		plotString <- paste0(plotString,")")
		eval(parse(text=plotString))
			##  changes to side panel
			##  changes to main panel
		if (input$plotChrIn == "") {
			insertUI(selector='#plotMain',where="afterBegin",ui=plotOutput('plotGenomeOut',width=990,height=1200))
		} else {
			insertUI(selector='#plotMain',where="afterBegin",ui=plotOutput('plotGenomeOut',width='70%',height=300))
		}
		insertUI(selector='#plotGenomeOut',where="afterEnd",ui=tags$div(id = 'downloadButtons'))
		insertUI(selector='#downloadButtons',ui=downloadButton("downloadPlot","Download plot"))
		insertUI(selector='#downloadButtons',ui=downloadButton("downloadCSV","Download data"))

	})

	observeEvent(input$smoothButton, {
		req(DFs$ratios,input$group,input$split)
		regionChr <- input$plotChrIn
		DFs$ratios<-smoothRatio(DFs$ratios,groupMin=as.integer(input$group),split=as.integer(input$split))
		updateSelectInput(session,'plotChrIn',selected=regionChr)
		for (i in 1:values$i) {
			removeUI(selector=paste0("#rmSmooth",i))
			insertUI(selector=paste0("#smoothie",i),where="afterBegin",ui=HTML(paste0(
				"<div id='smooth>",i,"'>Smooth</div>"
			)))
		}
	})

#########################################  outputs  #############################################

	output$plotGenomeOut <- renderPlot(plots$genomePlot)

	output$plotChrOut <- renderUI({
		chrList <- levels(DFs$ratios$chrom)
		selectInput("plotChrIn",NULL,c("Chr"="",chrList))
	})

	output$downloadCSV <- downloadHandler(
	    filename = function() {
			if (toggles$plotRegion==F) {
				paste0('ratiosData.csv')
			} else {
				region <- paste0(input$plotChrIn,"_",input$plotChrStart,"-",input$plotChrEnd)
				paste0('ratiosData_',region,'.csv')
			}
		},
	    content = function(file) {
			fileContent <- data.frame(
				chrom=factor(),chromStart=integer(),chromEnd=integer(),
				name.rep=factor(),name.nonRep=factor(),ratio=numeric(),
				group=integer(),splineSmooth=numeric()
			)
			guide <- DFs$guide
			guide <- guide[order(guide$order,na.last = T),]
			rownames(guide) <- 1:nrow(guide)
			samples <- length(na.omit(guide$order))
			for (i in 1:samples) {
				rep <- guide$name.rep[i]
				nonRep <- guide$name.nonRep[i]
				currentRatio <- DFs$ratios[DFs$ratios$name.rep==rep & DFs$ratios$name.nonRep==nonRep,]
				fileContent <- rbind(fileContent,currentRatio)
			warning(dfHead(fileContent))
			}
			if (toggles$plotRegion==F) {
				write.table(fileContent,file=file,sep=",",col.names=T,row.names=F,quote=F)
			} else {
				fileContent <- subset(
					fileContent,chrom==input$plotChrIn & chromStart>=input$plotChrStart & chromEnd<=input$plotChrEnd
				)
				write.table(fileContent,file=file,sep=",",col.names=T,row.names=F,quote=F)
			}
		}
	)

	output$downloadPlot <- downloadHandler(
		filename = function() {
			if (toggles$plotRegion==F) {
				paste0('ratiosPlot.pdf')
			} else {
				region <- paste0(input$plotChrIn,"_",input$plotChrStart,"-",input$plotChrEnd)
				paste0('ratiosPlot_',region,'.pdf')
			}
		},
		content = function(file) {
			if (toggles$plotRegion==F) {
				ggsave(file, plot=plots$genomePlot, device=cairo_pdf, width = 35, height = 50, units = "cm")
			} else {
				ggsave(file, plot=plots$genomePlot, device=cairo_pdf, width = 30, height = 12, units = "cm")
			}
		}
	)






#################################################~~~~~~~~~~~~~###############################################################
#################################################  STATS TAB  ###############################################################
#################################################~~~~~~~~~~~~~###############################################################


		##  STATS: load sample data
	observeEvent(input$exampleStats, {  ##  change relevant toggles and populate DFs$stats
		toggles$isExampleStats <- T
		DFs$stats <- sortSeq
	},ignoreInit=T)

		# add select ratio ui
	observeEvent(toggles$statSelectIn, {
			insertUI(selector='#statsSide',where="afterBegin",ui=HTML("<div id='statsSelectDiv'></div>"))
			insertUI(selector='#statsSelectDiv',where="afterBegin",ui=uiOutput("firstRatioOut"))
			insertUI(selector='#statsSelectDiv',where="beforeEnd",ui=uiOutput("secondRatioOut"))
	},ignoreInit=T)

		##  STATS: render select ratio ui
	output$firstRatioOut <- renderUI({
		selectInput('firstRatioIn',"Select ratios to calculate stats:",c("First ratio"="",names$ratio$name),width="75%")
	})

	output$secondRatioOut <- renderUI({
		selectInput("secondRatioIn",NULL,c("Second ratio"="",names$ratio$name),width="75%")
	})

	observeEvent(input$firstRatioIn, {
		if (input$firstRatioIn != "") {
			updateSelectInput(session,'secondRatioIn',choices=c("Second ratio"="",names$ratio$name[-which(names$ratio$name==input$firstRatioIn)]),selected=input$secondRatioIn)
		} else {
			updateSelectInput(session,'secondRatioIn',choices=c("Second ratio"="",names$ratio$name),selected=input$secondRatioIn)
		}
	})

	observeEvent(input$secondRatioIn, {
		if (input$secondRatioIn != "") {
			updateSelectInput(session,'firstRatioIn',choices=c("First ratio"="",names$ratio$name[-which(names$ratio$name==input$secondRatioIn)]),selected=input$firstRatioIn)
		} else {
			updateSelectInput(session,'firstRatioIn',choices=c("First ratio"="",names$ratio$name),selected=input$firstRatioIn)
		}
	})

	observeEvent({
		input$firstRatioIn
		input$secondRatioIn
	} ,{
		if (input$firstRatioIn != "" & input$secondRatioIn != "") {
			insertUI(selector='#statsSelectDiv',where="beforeEnd",ui=actionButton('calcStats',"Calculate"))
		} else {
			removeUI(selector='#calcStats')
		}
	})

	observeEvent(input$calcStats, {
		ratios <- subset(names$ratio,name==input$firstRatioIn)
		ratios[2,] <-  subset(names$ratio,name==input$secondRatioIn)
		DFs$stats <- ratioStats(
			subset(DFs$ratios,name.rep==ratios[1,1] & name.nonRep==ratios[1,2]),
			subset(DFs$ratios,name.rep==ratios[2,1] & name.nonRep==ratios[2,2]),
			names = ratios$name
		)
		removeUI('#calcStats')
	})

	observeEvent(DFs$stats, {
		ratioNames <- unique(DFs$stats[,c(as.character("name.rep"),as.character("name.nonRep"))])
		ratioNames$name <- paste0(ratioNames$name.rep," (",ratioNames$name.nonRep,")")
		if (toggles$statsCtrlsIn == F) {
				## changes to the STATS side panel
			removeUI(selector='#statsDescription')
			removeUI(selector='#exampleStatsDiv')
			insertUI(selector='#statSamples',where="afterBegin",ui=HTML("
				<div class='myTooltip' style='padding-bottom:15px;'><ctrlH>Samples</ctrlH>
					<span class='myTooltiptext'>This area contains individual sample controls</span>
				</div>
				<div id='statsSampleArea' style='padding-left:10px;'></div>
			"))
			insertUI(
				selector='#statSamples',
				where="afterEnd",
				ui=HTML("
					<div id='statsPlotElements'>
						<div class='myTooltip'><ctrlH>Additional features</ctrlH>
							<span class='myTooltiptext'>Enchance plot with your data (must be in bed format)</span>
						</div>
						<div style='padding:15px 10px;'>
							<div id='statSelectElement' class='inline' style='width:35%;'></div>
							<div id='statsUploadElement' class='inline' style='width:60%;padding-left:15px;'></div>
							<div id='statsPlotFeaturesDiv'></div>
						</div>
					</div>

					<div>
						<div class='myTooltip'><ctrlH>Plotting controls</ctrlH>
							<span class='myTooltiptext'>Choose the appearance of the plot</span>
						</div>
						<div id='statsPlottingCtrls' style='padding:10px;'>
							<div id='statsPlotRegionDiv' style='padding-bottom:10px;'>
								<div id='statsPlotRegion'>
									<b>Select region to plot</b><br>
									<div>
										<div id='statsPlotChrDiv' class='inline' style='width:20%;'></div>
										<div id='statsPlotChrStartDiv' class='inline' style='width:23%;'></div>
										<div id='statsPlotChrEndDiv' class='inline' style='width:23%;'></div>
										<div id='statsResetPlotRegionDiv' class='inline' style='width:28%;'></div>
									</div>
								</div>
							</div>
							<div id='statsGeomInputDiv' class='inline' style='width:45%;'>
								<div style='padding-top:5px;'><b>Plot type:</b></div>
								<div id='statsGeomDiv' style='padding-right:10px;width:80%;'></div>
							</div>
							<div id='statsPlotLimits' class='inline' style='width:45%;'>
								<div style='padding-top:5px;'><b>y limits:</b></div>
								<div>
									<div id='statsPlotLim1' class='inline' style='width:45%;'></div>
									<div id='statsPlotLim2' class='inline' style='width:45%;'></div>
								</div>
							</div>
							<div id='statsPlotGenomeBtnDiv' style='padding-top:20px;text-align:center;'></div>
						</div>
					</div>
				")
			)

			insertUI(selector="#statSelectElement",where="afterBegin",ui=selectInput('statsPlotFeatures',NULL,
				c("Choose type:"="","vLine"="Lines","Circle"="Circles","Rectangle"="Rectangles","Pointer"="Pointers")))
			insertUI(selector="#statsUploadElement",where="afterBegin",ui=fileInput(
				'statsPlotFeaturesFile',NULL,multiple=F,accept=".bed",buttonLabel = "Browse...", placeholder="No file selected"
			))

			insertUI(selector="#statsPlotChrDiv",where="afterBegin",ui=uiOutput('statsPlotChrOut',width="90%"))
			insertUI(
				selector="#statsPlotChrStartDiv",
				where="afterBegin",
				ui=textInput('statsPlotChrStart',NULL,NULL,width="95%",placeholder='Start')
			)
			insertUI(
				selector="#statsPlotChrEndDiv",
				where="afterBegin",
				ui=textInput('statsPlotChrEnd',NULL,NULL,width="95%",placeholder='End')
			)
			insertUI(selector='#statsResetPlotRegionDiv',where="afterBegin",ui=actionButton("statsResetPlotRegion","Reset region"))

			insertUI(selector="#statsGeomDiv",where="afterBegin",ui=selectInput(
				'statsGeomInput',
				NULL,
				c("scatter"='geom_point',"polygon"='geom_ribbon',"bar"='geom_segment'),
				selected='scatter')
			)
			insertUI(
				selector="#statsPlotLim1",
				where="afterBegin",
				ui=textInput('statsPlotLimLow',NULL,1,width="95%",placeholder='From')
			)
			insertUI(
				selector="#statsPlotLim2",
				where="afterBegin",
				ui=textInput('statsPlotLimHi',NULL,2,width="95%",placeholder='To')
			)
			insertUI(selector="#statsPlotGenomeBtnDiv",where="afterBegin",ui=actionButton('statsPlotGenomeButton'," Plot "))


				##  Populate the samples div
			insertUI(selector='#statsSampleArea',where="beforeEnd",ui=HTML(paste0("
				<div id='statsRatioDiv1' style='color:#696969;padding-bottom:10px;'>
					<div id='statsSampleName1' style='padding-bottom:10px;color:#000000;font-style:italic;font-weight: bold'>
						<div id='statsSample1'>",
							ratioNames$name.rep[1]," (",ratioNames$name.nonRep[1],")
						</div>
					</div>
					<div id='statsOrderDiv1' class='inline' style='width:15%;'> </div>
					<div class='inline'><b style='display:block;padding-left:10px;'>Data</b>
						<div class='inline'>
							<div id='statsRawDiv1' class='inline' style='padding-left:10px;'></div>
							<div class='inline' style='display:inline-block;padding-top:12px;'>Raw</div>
						</div>
						<div class='inline'>
							<div id='statSmoothDiv1' class='inline' style='padding-left:12px;'></div>
							<div id='statSmoothie1' class='inline' style='display:inline-block;padding-top:12px;'>
								<div id='statsRmSmooth1'>Smooth</div>
							</div>
						</div>
					</div>
					<div class='inline'><b style='display:block;padding-left:10px;'>Colour</b>
						<div id='statsColorDiv1' style='padding-top:5px;padding-left:10px;'> </div>
					</div>
				</div>
				<div id='statsRatioDiv2' style='color:#696969;padding-bottom:10px 0;'>
					<div id='statsSampleName2' style='padding-bottom:10px;color:#000000;font-style:italic;font-weight: bold'>
						<div id='statsSample2'>",
							ratioNames$name.rep[2]," (",ratioNames$name.nonRep[2],")
						</div>
					</div>
					<div id='statsOrderDiv2' class='inline' style='width:15%;'> </div>
					<div class='inline'><b style='display:block;padding-left:10px;'>Data</b>
						<div class='inline'>
							<div id='statsRawDiv2' class='inline' style='padding-left:10px;'></div>
							<div class='inline' style='display:inline-block;padding-top:12px;'>Raw</div>
						</div>
						<div class='inline'>
							<div id='statSmoothDiv2' class='inline' style='padding-left:12px;'></div>
							<div id='statSmoothie2' class='inline' style='display:inline-block;padding-top:12px;'>
								<div id='statsRmSmooth2'>Smooth</div>
							</div>
						</div>
					</div>
					<div class='inline'><b style='display:block;padding-left:10px;'>Colour</b>
						<div id='statsColorDiv2' style='padding-top:5px;padding-left:10px;'> </div>
					</div>
				</div>")
				)
			)
			insertUI(selector='#statsOrderDiv1',where="afterBegin",ui=selectInput('statsOrderInput1',"Order",c(1,2),selected=1))
			insertUI(selector='#statsOrderDiv2',where="afterBegin",ui=selectInput('statsOrderInput2',"Order",c(1,2),selected=2))
			insertUI(selector='#statsRawDiv1',where="afterBegin",ui=checkboxInput('statsRawInput1',NULL,TRUE))
			insertUI(selector='#statsRawDiv2',where="afterBegin",ui=checkboxInput('statsRawInput2',NULL,TRUE))
			insertUI(selector='#statSmoothDiv1',where="afterBegin",ui=checkboxInput('statSmoothInput1',NULL,FALSE))
			insertUI(selector='#statSmoothDiv2',where="afterBegin",ui=checkboxInput('statSmoothInput2',NULL,FALSE))
			insertUI(selector='#statsColorDiv1',where="afterBegin",ui=colourpicker::colourInput('statsColorInput1',NULL,"#7F7F7F"))
			insertUI(selector='#statsColorDiv2',where="afterBegin",ui=colourpicker::colourInput('statsColorInput2',NULL,"#00688B"))

				## changes to the STATS main panel
			insertUI(selector='#statsMain',where="afterBegin",ui=HTML("<div id='statsPlot'> <div>"))
			insertUI(selector='#statsPlot',where="afterBegin",ui=plotOutput('statsPlotOut',width=990,height=1200))
			insertUI(selector='#statsPlot',where="afterEnd",ui=HTML("<div id='statsDownloadButtons'></div>"))
			insertUI(selector='#statsDownloadButtons',ui=downloadButton("downloadStatsPlot","Download plot"))
			insertUI(selector='#statsDownloadButtons',ui=downloadButton("downloadStatsCSV","Download data"))
			toggles$statsCtrlsIn <- T
		} else {
			removeUI(selector='#statsSample1')
			removeUI(selector='#statsSample1')
			insertUI(
				selector='#statsSampleName1',
				where="afterBegin",
				ui=HTML(paste0("<div id='statsSample1'>",ratioNames$name.rep[1]," (",ratioNames$name.nonRep[1],")</div>")
			))
			insertUI(
				selector='#statsSampleName2',
				where="afterBegin",
				ui=HTML(paste0("<div id='statsSample2'>",ratioNames$name.rep[2]," (",ratioNames$name.nonRep[2],")</div>")
			))
		}
			##  construct guide DF
		DFs$guide = data.frame(
			order=as.integer(c(1,2)),
			name.rep=ratioNames$name.rep,
			name.nonRep=ratioNames$name.nonRep,
			raw=as.logical(c(TRUE,TRUE)),
			smooth=as.logical(c(FALSE,FALSE)),
			color=as.character(c("#7F7F7F","#00688B")),
			stringsAsFactors=F
		)
			##  Initial plot
		plots$statsPlot <- plotGenome(DFs$stats,guide=DFs$guide,plotting=F,geom="geom_point")
	},ignoreInit=T)

	output$statsPlotChrOut <- renderUI({
		chrList <- levels(DFs$stats$chrom)
		selectInput("statsPlotChrIn",NULL,c("Chr"="",chrList))
	})

	output$statsPlotOut <- renderPlot(plots$statsPlot)

	observeEvent(input$statsPlotFeaturesFile, {
		if (input$statsPlotFeatures!="") {
			filePath <- input$statsPlotFeaturesFile$datapath
			feature_str <- paste0("DFs$",input$statsPlotFeatures,"<-loadBed(filePath)")
			eval(parse(text=feature_str))
			str <- paste0("featureName <- levels(DFs$",input$statsPlotFeatures,"$name)[1]")
			eval(parse(text=str))

			features <- c("Lines","Circles","Rectangles","Pointers")
			colours <- c("#00FF00","#FFFFFF","#FF0000","#FF7F00")

			removeUI(paste0("#stats",input$statsPlotFeatures))

			insertUI('#statsPlotFeaturesDiv',"beforeEnd",ui=HTML(paste0("
				<div class='inline'>
					<div id='stats",input$statsPlotFeatures,"' class='inline' style='max-width:75%;padding-left:10px;'></div>
					<div id='statsColour",input$statsPlotFeatures,"' class='inline' style='max-width:25%;padding-left:10px;'></div>
				</div>"
			)))
			insertUI(paste0("#stats",input$statsPlotFeatures),"afterBegin",checkboxInput(paste0("stats",input$statsPlotFeatures),paste(featureName),TRUE))

			insertUI(paste0("#statsColour",input$statsPlotFeatures),"afterBegin",colourpicker::colourInput(
				paste0("stats",input$statsPlotFeatures,"ColourInput"),NULL,colours[match(input$statsPlotFeatures,features)],showColour = "background",palette = "limited"))
		}
	})

	observeEvent(input$statsPlotChrIn, {
		if (input$statsPlotChrIn!="") {
			chrom <- input$statsPlotChrIn
			chromEnd <- max(DFs$stats$chromEnd[DFs$stats$chrom == chrom])
			updateTextInput(session,'statsPlotChrStart',value=0)
			updateTextInput(session,'statsPlotChrEnd',value=chromEnd)
			toggles$statsRegion <- T
		}
	})

	observeEvent(input$statsResetPlotRegion, {
		updateSelectInput(session,'statsPlotChrIn',selected="")
		updateTextInput(session,'statsPlotChrStart',value="")
		updateTextInput(session,'statsPlotChrEnd',value="")
		toggles$statsRegion <- F
	})

	observeEvent(input$statsPlotGenomeButton, {
		req(DFs$stats,input$statsPlotLimLow,input$statsPlotLimHi)
		removeUI(selector='#statsPlotOut')
		ratioNames <- unique(DFs$stats[,c("name.rep","name.nonRep")])
		DFs$guide$name.rep <- as.character(ratioNames$name.rep)
		DFs$guide$name.nonRep <- as.character(ratioNames$name.nonRep)

	##  direct action: read inputs
		DFs$guide$order <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.integer(input[[paste0("statsOrderInput",i)]])
		})) #order
		DFs$guide$raw <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.logical(input[[paste0("statsRawInput",i)]])
		})) #raw
		DFs$guide$smooth <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.logical(input[[paste0("statSmoothInput",i)]])
		})) #smooth
		DFs$guide$color <- unlist(lapply(1:nrow(DFs$guide),function(i) {
			as.character(input[[paste0("statsColorInput",i)]])
		})) #color
		ylims <- c(as.numeric(input$statsPlotLimLow),as.numeric(input$statsPlotLimHi))

	##  direct action: create plot object
		plotString <- "plots$statsPlot <- plotGenome(DFs$stats,plotting=F,guide=DFs$guide,geom=input$statsGeomInput,ylims=ylims"
		if (input$statsPlotChrIn!="") {
			insertUI(selector='#statsPlot',where="afterBegin",ui=plotOutput('statsPlotOut',width='70%',height=300))
			plotString <- paste0(plotString,",region='",input$statsPlotChrIn,":",input$statsPlotChrStart,"-",input$statsPlotChrEnd,"'")
		} else {
			insertUI(selector='#statsPlot', where="afterBegin", ui=plotOutput('statsPlotOut',width=990,height=1200))
		}
		if (!is.null(DFs$Lines)) {
			if (input$statsLines) {
				plotString <- paste0(plotString,",lines=DFs$Lines,colourLines=input$statsLinesColourInput")
			}
		}
		if (!is.null(DFs$Circles)) {
			if (input$statsCircles) {
				plotString <- paste0(plotString,",circles=DFs$Circles,colourCircles=input$statsCirclesColourInput")
			}
		}
		if (!is.null(DFs$Rectangles)) {
			if (input$statsRectangles) {
				plotString <- paste0(plotString,",rectangles=DFs$Rectangles,colourRectangles=input$statsRectanglesColourInput")
			}
		}
		if (!is.null(DFs$Pointers)) {
			if (input$statsPointers) {
				plotString <- paste0(plotString,",pointers=DFs$Pointers,colourPointers=input$statsPointersColourInput")
			}
		}
		plotString <- paste0(plotString,")")
		eval(parse(text=plotString))

	##  changes to side panel
	##  changes to main panel
	})

	output$downloadStatsPlot <- downloadHandler(
		filename = function() {
			if (toggles$statsRegion==F) {
				paste0('statsPlot.pdf')
			} else {
				region <- paste0(input$statsPlotChrIn,"_",input$statsPlotChrStart,"-",input$statsPlotChrEnd)
				paste0('statsPlot_',region,'.pdf')
			}
		},
		content = function(file) {
			if (toggles$statsRegion==F) {
				ggsave(file, plot=plots$statsPlot, device=cairo_pdf, width = 35, height = 50, units = "cm")
			} else {
				ggsave(file, plot=plots$statsPlot, device=cairo_pdf, width = 30, height = 12, units = "cm")
			}
		}
	)

	output$downloadStatsCSV <- downloadHandler(
	    filename = function() {
			if (toggles$statsRegion==F) {
				paste0('statsData.csv')
			} else {
				region <- paste0(input$statsPlotChrIn,"_",input$statsPlotChrStart,"-",input$statsPlotChrEnd)
				paste0('statsData_',region,'.csv')
			}
		},
	    content = function(file) {
			fileContent <- data.frame(
				chrom=factor(),chromStart=integer(),chromEnd=integer(),
				name.rep=factor(),name.nonRep=factor(),ratio=numeric(),
				group=integer(),splineSmooth=numeric(),p.value=numeric()
			)
			guide <- DFs$guide
			guide <- guide[order(guide$order,na.last = T),]
			rownames(guide) <- 1:nrow(guide)
			samples <- length(na.omit(guide$order))
			for (i in 1:samples) {
				rep <- guide$name.rep[i]
				nonRep <- guide$name.nonRep[i]
				currentRatio <- DFs$stats[DFs$stats$name.rep==rep & DFs$stats$name.nonRep==nonRep,]
				fileContent <- rbind(fileContent,currentRatio)
			warning(dfHead(fileContent))
			}
			if (toggles$statsRegion==F) {
				write.table(fileContent,file=file,sep=",",col.names=T,row.names=F,quote=F)
			} else {
				fileContent <- subset(
					fileContent,chrom==input$statsPlotChrIn & chromStart>=input$statsPlotChrStart & chromEnd<=input$statsPlotChrEnd
				)
				write.table(fileContent,file=file,sep=",",col.names=T,row.names=F,quote=F)
			}
		}
	)
}
