
DATPATH=./Data
FUNPATH=./Programs
RESPATH=./Results

all:           data figures simulations

data:          $(DATPATH)/pbc_data.Rdata
               
simulations:   $(RESPATH)/results.Rdata \
               $(RESPATH)/beta_stacy_post.Rdata \
               $(RESPATH)/simulations.Rdata 
               
figures:       $(RESPATH)/Figure.pdf \
               $(RESPATH)/Supplementary_Fig_1.pdf \
               $(RESPATH)/Supplementary_Fig_2.pdf

clean:
	rm -fv $(DATPATH)/*.Rdata
	rm -fv $(RESPATH)/*.*
	rm -fv $(FUNPATH)/*.Rout
	rm -fv ./*.Rout
	rm -fv ./nohup.out

cleanlogs:
	rm -fv $(FUNPATH)/*.Rout
	rm -fv ./*.Rout
	rm -fv ./nohup.out


####### Datasets


$(DATPATH)/pbc_data.Rdata: $(FUNPATH)/prep_data.R
	R CMD BATCH $(FUNPATH)/prep_data.R


####### Simulations

$(DATPATH)/simulations.Rdata: $(FUNPATH)/BSB_code.R \
													 $(FUNPATH)/functions.R \
													 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/BSB_code.R

$(DATPATH)/beta_stacy_post.Rdata: $(FUNPATH)/beta_stacy_post.R \
													 $(FUNPATH)/functions.R \
													 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/beta_stacy_post.R

$(DATPATH)/results.Rdata:  $(FUNPATH)/results.R \
											     $(DATPATH)/beta_stacy_post.Rdata \
												   $(DATPATH)/simulations.Rdata \
												   $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/results.R

####### Figures

$(RESPATH)/Figure.pdf: $(FUNPATH)/make_plots.R \
											 $(FUNPATH)/functions.R \
											 $(DATPATH)/results.Rdata \
											 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/make_plots.R

####### Supplementary analyses and results

$(RESPATH)/Supplementary_Fig_1.pdf \
$(RESPATH)/Supplementary_Fig_2.pdf: $(FUNPATH)/supplementary_figures.R \
											 $(FUNPATH)/functions.R \
											 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/supplementary_figures.R