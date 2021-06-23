library(nearBynding)
# for supplemental figure, HNRNPK IP only, input only, and IP minus input for
# HepG2 and K562 cells

visualizeCapRStereogene(protein_file = c("ENCFF198ISB_liftOver",
                                         "ENCFF553XCL_liftOver"),
                        protein_file_input = "ENCFF019JFZ_liftOver",
                        CapR_prefix = "HepG2_3UTR",
                        legend = F,
                        out_file = "../HNRNPK_HepG2",
                        y_lim = c(-1.8, 3))
visualizeCapRStereogene(protein_file = c("ENCFF198ISB_liftOver",
                                         "ENCFF553XCL_liftOver"),
                        CapR_prefix = "HepG2_3UTR",
                        legend = F,
                        out_file = "../HNRNPK_HepG2_foreground",
                        y_lim = c(-1.8, 3))
visualizeCapRStereogene(protein_file = "ENCFF019JFZ_liftOver",
                        CapR_prefix = "HepG2_3UTR",
                        legend = F,
                        out_file = "../HNRNPK_HepG2_background",
                        y_lim = c(-1.8, 3))
visualizeCapRStereogene(protein_file = c("ENCFF405ESF_liftOver",
                                         "ENCFF894NKS_liftOver"),
                        protein_file_input = "ENCFF399CEH_liftOver",
                        CapR_prefix = "K562_3UTR",
                        legend = F,
                        out_file = "../HNRNPK_K562",
                        y_lim = c(-1.8, 3))
visualizeCapRStereogene(protein_file = c("ENCFF405ESF_liftOver",
                                         "ENCFF894NKS_liftOver"),
                        CapR_prefix = "K562_3UTR",
                        legend = F,
                        out_file = "../HNRNPK_K562_foreground",
                        y_lim = c(-1.8, 3))
visualizeCapRStereogene(protein_file = "ENCFF399CEH_liftOver",
                        CapR_prefix = "K562_3UTR",
                        legend = F,
                        out_file = "../HNRNPK_K562_background",
                        y_lim = c(-1.8, 3))
