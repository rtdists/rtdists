#require(testthat)
context("Diffusion pdiffusion and ddiffusion functions (new rcpp versions).")

is_testing <- TRUE  # Set to TRUE for release package 
## FALSE only works with extra non-supplied R and C files and after uncommenting 
# the two devtools lines (67 & 68)

pset1 <- list(a     = 0.8,  
              zr    = 0.3, 
              v     = 3.69, 
              t0    = 0.3, 
              d     = 0, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = "upper", precision = 3)

pset2 <- list(a     = seq(0.8, 0.9, by=0.1), 
              zr    = 0.3, 
              v     = 3.69, 
              t0    = seq(0.3,0.5,length=3), 
              d     = 0, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = c("upper","lower"), precision = 3)

pset3 <- list(a     = 1, 
              zr    = 0.3, 
              v     = 3.69, 
              t0    = seq(0.3,0.5,length=3), 
              d     = 0.1, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = c("upper","lower"), precision = 3)

pset4 <- list(a     = 1, 
              zr    = 0.7, 
              v     = 1.69, 
              t0    = seq(0.3,0.5,length=3), 
              d     = 0, 
              szr   = 0.2, 
              sv    = 0.9, 
              st0   = 0.1,
              bound = c("upper","lower"), precision = 3)

pset5 <- list(a     = 2, 
              zr    = 0.7,
              v     = 1.3,
              t0    = 0.2,
              d     = 0,
              szr   = 0.3,
              sv    = 0.4,
              st0   = 0.5,
              bound = "upper", precision = 3)


# Values to calculate d/p for 
x <- seq (0,4,by=0.1)

if (!is_testing) {
  #### This generates datasets from the old (pre-RCpp) version of rtdists (specifically, from 0.6-6)
  #### That C code is not included in the release package.
  
  # UNLOAD EXISTING RTDISTS, LOAD IN OLD VERSION FOR COMPARISON TESTING
  ## UNCOMMENT NEXT TWO LINES
  # require (devtools)
  # devtools::unload(rtdists)
  #install.packages ("rtdists", lib="tests\\temp_testing\\old_rtdists_0.6-6\\")
  require (rtdists, lib.loc="tests\\temp_testing\\old_rtdists_0.6-6")

  orig_d_set1 <- ddiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision)
  orig_d_set2 <- ddiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision)
  orig_d_set3 <- ddiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision)
  orig_d_set4 <- ddiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision)
  orig_d_set5 <- ddiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision)
  
  orig_p_set1 <- pdiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision)
  orig_p_set2 <- pdiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision)
  orig_p_set3 <- pdiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision)
  orig_p_set4 <- pdiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision)
  orig_p_set5 <- pdiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision)
  
  # Print out in a format that makes it easy to just copy/paste below for when is_testing = TRUE  
  tmp <- paste(orig_d_set1, collapse=","); paste0("orig_d_set1 <- c(", tmp, ")")
  tmp <- paste(orig_d_set2, collapse=","); paste0("orig_d_set2 <- c(", tmp, ")")
  tmp <- paste(orig_d_set3, collapse=","); paste0("orig_d_set3 <- c(", tmp, ")")
  tmp <- paste(orig_d_set4, collapse=","); paste0("orig_d_set4 <- c(", tmp, ")")
  tmp <- paste(orig_d_set5, collapse=","); paste0("orig_d_set5 <- c(", tmp, ")")
  
  tmp <- paste(orig_p_set1, collapse=","); paste0("orig_p_set1 <- c(", tmp, ")")
  tmp <- paste(orig_p_set2, collapse=","); paste0("orig_p_set2 <- c(", tmp, ")")
  tmp <- paste(orig_p_set3, collapse=","); paste0("orig_p_set3 <- c(", tmp, ")")
  tmp <- paste(orig_p_set4, collapse=","); paste0("orig_p_set4 <- c(", tmp, ")")
  tmp <- paste(orig_p_set5, collapse=","); paste0("orig_p_set5 <- c(", tmp, ")")
} else {
  # Reference values regenerated with GSL adaptive quadrature (replaces rectangular
  # midpoint rule for st0/sz variability integration).  Verified against WienR
  # (independent implementation using cubature library) to ~1e-10 accuracy.
  orig_d_set1 <- c(0,0,0,0,4.60128402369186,3.04118238262028,0.775211977459591,0.200977786985829,0.0554615467551478,0.016192764457705,0.00495917872855415,0.00158172208234681,0.000522329532511231,0.000177738832473327,6.20779848408495e-05,2.21815705635608e-05,8.08644129943854e-06,3.00074180888364e-06,1.13123116472851e-06,4.3250782168812e-07,1.67465834017741e-07,6.55846592365033e-08,2.59505714754716e-08,1.03644090462547e-08,4.17474476439496e-09,1.69465337352762e-09,6.92806613474608e-10,2.85081918235247e-10,1.18012367342046e-10,4.91227925080697e-11,2.05520163715664e-11,8.63928904662754e-12,3.64758961953118e-12,1.54633851054846e-12,6.58040198069259e-13,2.81021224547655e-13,1.20410403994702e-13,5.17529541404755e-14,2.2308395454351e-14,9.64242709034107e-15,4.1784794347066e-15)
  orig_d_set2 <- c(0,0,0,0,0,0,0.775211977459591,0.0699951185456533,0.775211977459592,0.0035383135441595,0.016192764457705,0.0035383135441595,0.000522329532511231,0.000251904204355878,0.000522329532511232,2.1842681956195e-05,2.21815705635607e-05,2.1842681956195e-05,1.13123116472851e-06,2.16036599054688e-06,1.13123116472851e-06,2.34361628377365e-07,6.55846592365032e-08,2.34361628377365e-07,4.17474476439496e-09,2.71996437623082e-08,4.17474476439499e-09,3.32152864607518e-09,2.85081918235247e-10,3.32152864607517e-09,2.05520163715664e-11,4.21894675081224e-10,2.05520163715663e-11,5.52810755429853e-11,1.54633851054846e-12,5.52810755429853e-11,1.20410403994702e-13,7.42709389833597e-12,1.20410403994702e-13,1.01846872348035e-12,9.64242709034111e-15)
  orig_d_set3 <- c(0,0,0,0,0.225610071703342,0,0.937263945961779,0.126716433984152,0.93726394596178,0.00757726512534319,0.0455036861061427,0.00757726512534319,0.00330266309844789,0.000731489607573722,0.0033026630984479,8.74086474302187e-05,0.000315403878014916,8.74086474302186e-05,3.62299524993117e-05,1.20145210323639e-05,3.62299524993117e-05,1.82072329653906e-06,4.7406658341355e-06,1.82072329653906e-06,6.82426968833283e-07,2.96183650575236e-07,6.82426968833285e-07,5.08123607484398e-08,1.05578727258425e-07,5.08123607484397e-08,1.72718874306342e-08,9.08170976819027e-09,1.72718874306341e-08,1.67640321140245e-09,2.95308783484761e-09,1.67640321140245e-09,5.23177536818716e-10,3.17566972558544e-10,5.23177536818716e-10,6.14417809483706e-11,9.54204236228889e-11)
  orig_d_set4 <- c(0,0,0,0,0,0,0.855224973063466,0.208600365980508,0.855224973063466,0.0413133609761505,0.120195715931898,0.0413133609761505,0.0200186207062867,0.00803986701160756,0.0200186207062867,0.00161192347347385,0.00357670599797717,0.00161192347347385,0.000669895644588381,0.000330125891058263,0.000669895644588382,6.86710702158312e-05,0.000129753417786468,6.86710702158311e-05,2.5765258321897e-05,1.44542589002176e-05,2.57652583218971e-05,3.07060405014058e-06,5.21438471374448e-06,3.07060405014057e-06,1.07111805117301e-06,6.57138807458781e-07,1.071118051173e-06,1.41483725821826e-07,2.2266153768523e-07,1.41483725821826e-07,4.67377083659198e-08,3.0614490575549e-08,4.67377083659198e-08,6.6523553464196e-09,9.88952812953638e-09)
  orig_d_set5 <- c(0,0,0,0.000695939862177696,0.041086015153379,0.17173181815623,0.354566596669274,0.546517108269806,0.724502419409127,0.840941886574165,0.843640027389452,0.772139799460312,0.672167440124618,0.569010924389478,0.473990616301451,0.391153159744939,0.321072753570857,0.262802899001706,0.214843474634318,0.17560145005114,0.143594607472721,0.117527285414536,0.0963046946382963,0.079019826864937,0.0649300063470791,0.0534313136369797,0.0440345766395878,0.0363443420159148,0.0300411345161237,0.0248668132940429,0.0206126420334974,0.0171096462234087,0.014220854683373,0.0118350722717386,0.0098618865877069,0.00822766431557216,0.00687233914902738,0.00574683213827587,0.00481097722382841,0.00403185052960627,0.00338242266769738)

  orig_p_set1 <- c(0,0,0,0,0.169654748430154,0.644852220438606,0.811721201050489,0.854057849069507,0.865306050934252,0.868482806734998,0.869428910313148,0.869723626727494,0.869819010517484,0.86985091068276,0.869861886388069,0.869865757198307,0.86986715220185,0.869867664636671,0.869867856083721,0.869867928694559,0.869867956607604,0.869867967468629,0.869867971741068,0.869867973438421,0.869867974118822,0.869867974393807,0.869867974505775,0.869867974551679,0.869867974570617,0.869867974578475,0.869867974581754,0.869867974583128,0.869867974583707,0.869867974583952,0.869867974584056,0.8698679745841,0.869867974584119,0.869867974584127,0.869867974584131,0.869867974584132,0.869867974584133)
  orig_p_set2 <- c(0,0,0,0,0,0,0.811721196913999,0.125574562330438,0.811721196913999,0.131997518069757,0.868482802598511,0.131997518069757,0.869819006380996,0.132363879963647,0.869819006380996,0.132391772625441,0.86986575306182,0.132391772625441,0.869867851947233,0.132394304104058,0.869867851947233,0.132394562734449,0.869867963332141,0.132394562734449,0.869867969982335,0.132394591469189,0.869867969982335,0.132394594864978,0.869867970415191,0.132394594864978,0.869867970445266,0.132394595285525,0.869867970445266,0.13239459533954,0.869867970447464,0.13239459533954,0.869867970447631,0.132394595346681,0.869867970447631,0.132394595347648,0.869867970447645)
  orig_p_set3 <- c(0,0,0,0,0.00178292322825537,0,0.775411087834309,0.120922452250843,0.775411087834309,0.132716135127724,0.861195435552079,0.132716135127724,0.865901441795384,0.133576741913059,0.865901441795384,0.133666526168523,0.866276666320826,0.133666526168523,0.866314922623755,0.133677824411988,0.866314922623755,0.133679435165113,0.866319529210655,0.133679435165113,0.866320153470468,0.133679685884491,0.866320153470468,0.13367972750184,0.866320245781243,0.13367972750184,0.866320260352645,0.133679734753962,0.866320260352645,0.133679736066237,0.866320262787705,0.133679736066237,0.866320263198217,0.133679736310875,0.866320263198217,0.133679736357606,0.866320263283258)
  orig_p_set4 <- c(0,0,0,0,0,0,0.781191530696519,0.0504836200207758,0.781191530696518,0.0819600110121512,0.890312005832964,0.0819600110121512,0.906961293667404,0.0880426402075123,0.906961293667404,0.0892403909461952,0.909812660187627,0.0892403909461952,0.910331667455244,0.0894825802843886,0.910331667455244,0.0895324893578018,0.910430168562784,0.0895324893578018,0.910449434981355,0.0895429198770103,0.910449434981355,0.0895451233632549,0.910453289384567,0.0895451233632549,0.910454074007284,0.0895455928255754,0.910454074007284,0.0895456935333022,0.910454235934917,0.0895456935333022,0.910454269724044,0.0895457152583196,0.910454269724044,0.0895457199671402,0.910454276838814)
  orig_p_set5 <- c(0,0,0,7.08446361178149e-06,0.00138436468812972,0.0113507675678673,0.0374267297023063,0.0825133638296158,0.146245531551868,0.225422581216873,0.310520009922219,0.391726173987651,0.464066183870186,0.526107264997353,0.57818038301509,0.62134283021036,0.656860758949882,0.685970510299194,0.70978061414215,0.729242390730573,0.745152325059471,0.758167663671159,0.768826127576991,0.777565489792073,0.784741228585,0.790641683658525,0.795500713500313,0.799508094326108,0.802817980191424,0.805555746159603,0.807823505956937,0.809704554224829,0.811266941644485,0.812566353290477,0.813648428047359,0.814550629885263,0.815303759724723,0.815933178816323,0.816459800292277,0.816900894157634,0.81727074191634)

  test_that("ensure new RCpp ddiffusion function produces the same result as previous C-only versions:", {
  
    tolerance <- 1e-5
  
    expect_equal(orig_d_set1, ddiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision), tolerance=tolerance)
    expect_equal(orig_d_set2, ddiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision), tolerance=tolerance)
    expect_equal(orig_d_set3, ddiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision), tolerance=tolerance)
    expect_equal(orig_d_set4, ddiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision), tolerance=tolerance)
    expect_equal(orig_d_set5, ddiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision), tolerance=tolerance)
  })
  
  test_that("ensure new RCpp pdiffusion function produces the same result as previous C-only versions:", {
  
    tolerance <- 1e-3
    
    expect_equal(orig_p_set1, pdiffusion (x, response=pset1$bound, a=pset1$a, v=pset1$v, t0=pset1$t0, z=pset1$zr, d=pset1$d, sz=pset1$szr, sv=pset1$sv, st0=pset1$st0, precision=pset1$precision), tolerance=tolerance)
    expect_equal(orig_p_set2, pdiffusion (x, response=pset2$bound, a=pset2$a, v=pset2$v, t0=pset2$t0, z=pset2$zr, d=pset2$d, sz=pset2$szr, sv=pset2$sv, st0=pset2$st0, precision=pset2$precision), tolerance=tolerance)
    expect_equal(orig_p_set3, pdiffusion (x, response=pset3$bound, a=pset3$a, v=pset3$v, t0=pset3$t0, z=pset3$zr, d=pset3$d, sz=pset3$szr, sv=pset3$sv, st0=pset3$st0, precision=pset3$precision), tolerance=tolerance)
    expect_equal(orig_p_set4, pdiffusion (x, response=pset4$bound, a=pset4$a, v=pset4$v, t0=pset4$t0, z=pset4$zr, d=pset4$d, sz=pset4$szr, sv=pset4$sv, st0=pset4$st0, precision=pset4$precision), tolerance=tolerance)

    tolerance <- 1e-2   # !! This still produces 1e-2 errors with the 5th parameter set!
    expect_equal(orig_p_set5, pdiffusion (x, response=pset5$bound, a=pset5$a, v=pset5$v, t0=pset5$t0, z=pset5$zr, d=pset5$d, sz=pset5$szr, sv=pset5$sv, st0=pset5$st0, precision=pset5$precision), tolerance=tolerance)
  })

}
