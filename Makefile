
DATPATH=./Data
FUNPATH=./Programs
RESPATH=./Results

all:           data figures simulations

data:          $(DATPATH)/pbc_data.Rdata
               
simulations:   $(RESPATH)/results.Rdata \
               $(RESPATH)/beta_stacy_post.Rdata \
               $(RESPATH)/simulations.Rdata \
               $(RESPATH)/sims_censoring_final.rds 
               
figures:       $(RESPATH)/Figure.pdf \
               $(RESPATH)/Supplementary_Fig_1.pdf \
               $(RESPATH)/Supplementary_Fig_2.pdf \
               $(RESPATH)/GvdVa.pdf \
               $(RESPATH)/Censoring_surv.pdf \
               $(RESPATH)/Censoring_rmst.pdf

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

$(RESPATH)/simulations.Rdata: $(FUNPATH)/BSB_code.R \
													 $(FUNPATH)/functions.R \
													 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/BSB_code.R

$(RESPATH)/beta_stacy_post.Rdata: $(FUNPATH)/beta_stacy_post.R \
													 $(FUNPATH)/functions.R \
													 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/beta_stacy_post.R

$(RESPATH)/results.Rdata:  $(FUNPATH)/results.R \
											     $(RESPATH)/beta_stacy_post.Rdata \
												   $(RESPATH)/simulations.Rdata \
												   $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/results.R

####### Figures

$(RESPATH)/Figure.pdf: $(FUNPATH)/make_plots.R \
											 $(FUNPATH)/functions.R \
											 $(RESPATH)/results.Rdata \
											 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/make_plots.R

####### Supplementary analyses and results

$(RESPATH)/Supplementary_Fig_1.pdf \
$(RESPATH)/Supplementary_Fig_2.pdf: $(FUNPATH)/supplementary_figures.R \
											 $(FUNPATH)/functions.R \
											 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/supplementary_figures.R
	
$(RESPATH)/Censoring_surv.pdf \
$(RESPATH)/Censoring_rmst.pdf: $(FUNPATH)/plots_censoring.R \
															 $(FUNPATH)/functions.R \
															 $(RESPATH)/sims_censoring_final.rds
	R CMD BATCH $(FUNPATH)/plots_censoring.R

$(RESPATH)/sims_censoring_final.rds: $(FUNPATH)/Censoring_rate_simulations.R \
																		 $(FUNPATH)/functions.R \
																		 $(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/Censoring_rate_simulations.R

$(RESPATH)/GvdVa.pdf: $(FUNPATH)/GvdVa_simulations.R \
											$(FUNPATH)/functions.R \
											$(DATPATH)/pbc_data.Rdata
	R CMD BATCH $(FUNPATH)/GvdVa_simulations.R