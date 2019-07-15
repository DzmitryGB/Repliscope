library(shiny)
library(ggplot2)
fluidPage(theme="Repliscope.css",
	headerPanel("Repliscope", windowTitle="Repliscope"),

sidebarLayout(
	sidebarPanel(
			conditionalPanel(condition="input.tabs=='About'",
				div(id="aboutSide",
					HTML("<h4 style='padding-bottom:10px;'><b>Selected publications</b></h4><ul style='text-align:justify;'>
							<li>
								Hawkins, M., Malla, S., Blythe, M. J., Nieduszynski, C. A., & Allers, T. (2013).
								<i>Accelerated growth in the absence of DNA replication origins.</i> Nature, 503(7477), 544–547.
								<a href='http://doi.org/10.1038/nature12650' target='_blank'>View</a>&nbsp;&nbsp;&nbsp;
								<a href='https://www.ncbi.nlm.nih.gov/pubmed/24185008' target='_blank'>Pubmed</a>
							</li>
							<br>
							<li>
								Natsume T, Müller CA, Katou Y, Retkute R, Gierliński M, Araki H, Blow JJ, Shirahige K, Nieduszynski CA, Tanaka TU. (2013).
								<i>Kinetochores coordinate pericentromeric cohesion and early DNA replication by Cdc7-Dbf4 kinase recruitment.</i> Mol Cell, 50(5):661-74.
								<a href='http://doi.org/10.1016/j.molcel.2013.05.011' target='_blank'>View</a>&nbsp;&nbsp;&nbsp;
								<a href='https://www.ncbi.nlm.nih.gov/pubmed/    23746350' target='_blank'>Pubmed</a>
							</li>
							<li>
								Müller,C.A., Hawkins,M., Retkute,R., Malla,S., Wilson,R., Blythe,M.J., Nakato,R., Komata,M., Shirahige,K.,
								de Moura,A.P.S., Nieduszynski, C.A. (2014). <i>The dynamics of genome replication using deep sequencing.</i>
								Nucleic Acids Res., 42, e3..&nbsp;&nbsp;&nbsp;
								<a href='http://doi.org/10.1093/nar/gkt878' target='_blank'>View</a>&nbsp;&nbsp;&nbsp;
								<a href='https://www.ncbi.nlm.nih.gov/pubmed/24089142' target='_blank'>Pubmed</a>
							</li>
							<br>
						</ul>
					")
				)
			),

	##~~~~~~~~~~~~~~~~~~~~~~~~~~  COVERAGE  ~~~~~~~~~~~~~~~~~~~~~~~~~~##


			conditionalPanel(condition="input.tabs=='Coverage'",
				div(id='coverageSide',
					div(id='loadBedDiv',
						class="inline",
						style="padding-bottom:30px;width:75%;margin-right:10px;",
						HTML("<div class='myTooltip'><label>Load a bed file:</label><span class='myTooltiptext'>
							File must be produced with the CANmapper script.
							</span></div>"),
						div(class="inline",
							fileInput("bedFile", NULL,multiple=F,accept=".bed",buttonLabel = "Browse...", placeholder="No file selected")
						)
					),
					div(id='analyseOrReset',
						class="inline",
						style="padding-top:25px;"
					),
					div(id='exampleCoverageDiv',
						style="margin-top:-10px;",
						HTML("<div style='padding-bottom:5px;'><b>Or</b></div>"),
						actionButton('exampleCoverage',"Load example")
					)
				)
			),

	##~~~~~~~~~~~~~~~~~~~~~~~~~~~  RATIO  ~~~~~~~~~~~~~~~~~~~~~~~~~~~##

			conditionalPanel(condition="input.tabs=='Ratio'",
				div(id='ratioSide',
					div(id='loadRatioDiv',
						style="margin-bottom:30px;width:75%;margin-right:10px;",
						HTML("<div class='myTooltip'><label>Load a saved ratio file:</label><span class='myTooltiptext'>
							File must be produced using this page
							</span></div>"),
						div(class='inline',
							fileInput("ratioFile",NULL,multiple=F,buttonLabel = "Browse...", placeholder="No file selected")
						)
					),
					div(id='exampleRatioDiv',
						style="margin-top:-10px;",
						HTML("<div style='padding-bottom:5px;'><b>Or</b></div>"),
						actionButton('exampleRatio',"Load example")
					)
				)
			),

	##~~~~~~~~~~~~~~~~~~~~~~~~~~~  PLOT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~##

			conditionalPanel(condition="input.tabs=='Plot'",
				div(id='plotSide',
					div(
						id='loadPlotDiv',
						style="margin-bottom:30px;width:75%;margin-right:10px;",
						HTML("<div class='myTooltip'><label>Load a saved ratios file:</label><span class='myTooltiptext'>
							File must be produced using this page</span></div>"
						),
						fileInput('plotFile',NULL,multiple=F,buttonLabel = "Browse...", placeholder="No file selected")
					),
					div(
						id='examplePlotDiv',
						style="margin-top:-10px;",
						HTML("<div style='padding-bottom:5px;'><b>Or</b></div>"),
						actionButton('examplePlot',"Load example")
					),
					div(
						id='plotSideCtrls',
						div(id='samples',style='padding-bottom:15px;')
					)
				)
			),



	##~~~~~~~~~~~~~~~~~~~~~~~~~~~  STATS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~##

			conditionalPanel(condition="input.tabs=='Stats'",
				div(id="statsSide",
					div(
						id='exampleStatsDiv',
						style='padding-top:30px;',
						actionButton('exampleStats',"Load example")
					),
					div(id='statSamples',style='padding-bottom:15px;')
				)
			),
			width=3
		),


		mainPanel(
			tabsetPanel(
				tabPanel("About",
					div(id="aboutMain",
						style="padding:20px;",
						HTML('
						<div class="col-sm-6"><h3>Replication time profiling using DNA relative copy number</h3>
						<p style="padding-top:10px;text-align:justify;color:#404040;">
							Repliscope is an R package for creating, normalising, comparing and plotting DNA replication timing profiles. The analysis pipeline
							starts with BED-formatted read count files (output of <a href="https://github.com/DzmitryGB/localMapper" target="_blank">localMapper</a>) obtained by
							high-throughput sequencing of DNA from replicating and non-replicating cells. There are three methods of measuring DNA replication
							dynamics using relative copy number (Fig): sort-seq, sync-seq and marker frequency analysis (MFA-seq). Sort-seq uses fluorescence-activated
							cell sorting (FACS) to enrich for non-replicating and replicating cells from an asynchronous population. Sync-seq requires cells to be
							arrested in non-replicating cell cycle phase (i.e. G1), followed by release into S phase. Samples are then taken throughout S phase
							when cells synchronously synthesise DNA according to the replication timing programme. In the case of MFA-seq, rapidly dividing cells in
							exponential growth phase are directly used as the replicating sample, while a saturated culture serves as a non-replicating control
							sample. While the latter approach of obtaining cells is the simplest, it also requires deeper sequencing due to decreased dynamic range
							and, thus, is more suitable for organisms with small genomes (typically, bacteria).
						</p></div>
						<div class="col-sm-6">
							<div style="overflow:hidden;max-width:767px;min-width:400px;padding:50px 10px;"><img src="outline.png" width="100%"></div>
						</div>')
					)
				),
				tabPanel("Coverage",
					div(id="coverageMain",
						HTML("
							<div id='coverageDescription' style='padding:50px;text-align:justify;color:#404040;'>
								<p>Use the menu on the left to either load example data or upload your own. The uploaded bed file may or may not have a
								header and <b>must</b> contain 5 columns (<b><i>chrom, chromStart, chromEnd, name, score</i></b>) where the first three columns
								define the genomic bin and the \"score\" column is used for storing coverage data (reads per bin). The \"name\"
								column will be repurposed to store the bed file's name, therefore bed file names should be descriptive.</p>
								<p>When the file is uploaded, a snippet of its modified content is displayed initially - please make sure everything looks right
								before continuing.</p>
								<p><b>Tip:</b> Hover over headings in the control panel to discover tooltips!</p>
							</div>
						")
					)
				),
				tabPanel("Ratio",
					div(id="ratioMain",
						HTML("
							<div id='ratioDescription' style='padding:50px;text-align:justify;color:#404040;'>
								<p>Use the menu on the left to either load example data or upload your own (created earlier using this page!).
								Once you have saved at least one sample for each replicating and non-replicating sample type, this page
								will display elements for making a new ratio. Initially, the ratio is normalised by the total read
								number and will have a distribution around one.</p>
                <p><b>Trimming</b> should be done if there are ratio values far outside of the main population, as they will skew the automatic
								normalisation. A range of 0.5-1.5 is a very safe starting point.</p>
								<p>In the case of full range S phase samples (sorted whole S phase or synchronised
								S phase population, where at least some regions are completely replicated), <b>automatic normalisation</b> may be used.
								This scales the data to lie between one and two, based on minimising the sum of data points outside of this region. Same strategy can be used
                with an asynchronous cell culture (marker frequency analysis), but the upper limit of the scale should be adjusted accordingly. For example, if the
                asynchronous population contains 20% of cells in S phase, the upper limit should be set to 1.2.</p>
								<p>If replicating samples come from S phase timepoints of a cell cycle experiment, <b>manual normalisation</b> should be used. The median values of DNA
                fluorescence obtained using flow cytometry of the samples during cell cycle experiment should be fitted to a sigmoidal function. We provide an <a
                href='https://dzmitry.shinyapps.io/flowfit/' target='_blank'>online tool</a> to simplify the process. The normalisation factor for the synchronised
                samples can be calculated as 1 + [median of replicated bulk DNA]. For example, if
                cells from a synchronised S phase population are have 15% bulk genome replication completed, a normalisation  factor of 1.15 should be used.</p>
							</div>
						")
					)
				),
				tabPanel("Plot",
					div(id="plotMain",
						HTML("
							<div id='plotDescription' style='padding:50px;text-align:justify;color:#404040;'>

								<p>Use menu on the left to either load an example data or upload your own multiple ratios file (created earlier
								using this page!). Any ratios saved in the previous tab should appear on the left, along with some controls.</p>

								<p><b>Additional features</b> may be added to the plot. Typically, replication origins are plotted as circles and centromere
								location as vertical lines. Rectangles can be used to highlight a big region of a chromosome, while pointers can be used to
								pinpoint a specific locus. Use the \"name\" field in the bed file to name the feature.</p>

								<p> Run <b>smoothing</b> to plot smoothed data. Due to the discrete nature of the data (separate chromosomes,
								missing bins), smoothing is done in groups. <label>Group size</label> controls a minimum number of datapoints
								to smooth (each group must contain a minimum of 4 datapoints for the spline algorithm to work).
								<label>Split</label> value controls the number of missing bins along a chromosome to initiate a new
								group.</p>

								<p>Use <b>plotting controls</b> to zoom into a particular chromosome (or a subchromosomal region) or
								choose different plot type (note that the only suitable plot type to display both raw and smooth data is scatter plot).
								Y axis limits allow to plot marker frequency analysis data.</p>
							</div>
						")
					)
				),
				tabPanel("Stats",
					div(id="statsMain",
						HTML("
							<div id='statsDescription' style='padding:50px;text-align:justify;color:#404040;'>
								<p>Use menu on the left to load example data.</p>
								<p>Any ratios saved in the 'Ratio' tab should appear on the left (when at least two are available). Plotting controls are
								similar to the 'Plot' tab.</p>
								<p>The statistical analysis is based on <a href='https://en.wikipedia.org/wiki/Standard_score' target='_blank'>z-scores</a>.
								Statistically significant differences between two replication profiles are separated into bins of 0.99-0.999 significance
								(p-value**) and above 0.999 significance (p-value***).</p>
							</div>
						")
					)
				),id="tabs"
			), width=9
		)
	),
	HTML("<div class='footer'><span style='height:50px;line-height:50px;vertical-align:middle;'>
		<a href='http://nieduszynski.org/' target='_blank'>Nieduszynski lab</a></span></div>")
)
