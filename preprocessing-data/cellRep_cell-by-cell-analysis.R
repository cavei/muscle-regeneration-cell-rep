
# Load the library.
library(cellCB)
library(RUVSeq)

# We create a directory to store the analysis.

resDir <- "rpkm-RDatas-padj0.01-lfc2/"
dir.create(resDir, showWarnings = F)

ec_wt <- runWithParametersRpkm("../expression-180520/countData.txt",
                               "../expression-180520/colData.txt",
                               "ec", "wt",
                               runWithK = 2, runWithFormula = "~days + W_1 + W_2",
                               pADJ.thr = 0.01, logfcthr = 2)

save(ec_wt, file=paste0(resDir,"ec","-","wt",".RData"))


fap_wt <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                                "fap", "wt", runWithK = 1, runWithFormula = "~days + W_1",
                                pADJ.thr = 0.01, logfcthr = 2)

save(fap_wt, file=paste0(resDir,"fap","-","wt",".RData"))

mp_wt <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                               "mp", "wt", runWithK = 1, runWithFormula = "~days + W_1",
                               pADJ.thr = 0.01, logfcthr = 2)
save(mp_wt, file=paste0(resDir,"mp","-","wt",".RData"))

inf_wt <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                                "inf", "wt", runWithK = 1, runWithFormula = "~days + W_1",
                                pADJ.thr = 0.01, logfcthr = 2)
save(inf_wt, file=paste0(resDir,"inf","-","wt",".RData"))

per_wt <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                                "per", "wt", runWithK = 1, runWithFormula = "~days + W_1",
                                removeSamples = "r1_3_per-wt", pADJ.thr = 0.01, logfcthr = 2)
save(per_wt, file=paste0(resDir,"per","-","wt",".RData"))

tot_wt <- runWithParametersRpkm("../full-muscle-180520/total-countData.txt",
                                     "../full-muscle-180520/t-colData.txt",
                                     "tot", "wt", runWithK = 1, runWithFormula = "~days + W_1",
                                     pADJ.thr = 0.01, logfcthr = 2)

save(tot_wt, file=paste0(resDir,"tot","-","wt",".RData"))

# We move to the KO condition.

ec_ko <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                               "ec", "ko", runWithK = 3, runWithFormula = "~days + W_1 + W_2 + W_3",
                               removeSamples = "r1_7_ec-ko",
                               pADJ.thr = 0.01, logfcthr = 2)
save(ec_ko, file=paste0(resDir,"ec","-","ko",".RData"))


fap_ko <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                                "fap", "ko", runWithK = 1, runWithFormula = "~days + W_1",
                                pADJ.thr = 0.01, logfcthr = 2)
save(fap_ko , file=paste0(resDir,"fap","-","ko",".RData"))

mp_ko <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                               "mp", "ko", runWithK = 1, runWithFormula = "~days + W_1",
                               removeSamples = "r1_14_mp-ko", pADJ.thr = 0.01, logfcthr = 2)
save(mp_ko, file=paste0(resDir,"mp","-","ko",".RData"))


all_wt <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                                "ALL", "wt", runWithK=1, runWithFormula="~cell + W_1", 
                                removeSamples = "r1_3_per-wt",
                                varName = "cell", pADJ.thr = 0.01, logfcthr = 2)
save(all_wt, file=paste0(resDir,"all","-","wt",".RData"))


all_ko <- runWithParametersRpkm("../expression-180520/countData.txt",
                                "../expression-180520/colData.txt",
                                "ALL", "ko", runWithK=1, runWithFormula="~cell + W_1",
                                removeSamples = c("r1_14_mp-ko","r1_7_ec-ko"), varName = "cell", pADJ.thr = 0.01, logfcthr = 2)
save(all_ko, file=paste0(resDir,"all","-","ko",".RData"))

all_wt_ko <- runWithParametersRpkm("../expression-180520/countData.txt",
                                   "../expression-180520/colData_plus.txt",
                                   "ALL", c("wt","ko"), runWithK=1, runWithFormula="~cell_condition + W_1",
                                   removeSamples = c("r1_3_per-wt", "r1_14_mp-ko", "r1_7_ec-ko"),
                                   varName = "cell_condition", pADJ.thr = 0.01, logfcthr = 2)

save(all_wt_ko, file=paste0(resDir,"all","-","wt", "-","ko",".RData"))
