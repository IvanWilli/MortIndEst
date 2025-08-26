# get a model lt from adjusted qx with Brass method using census CEB and CS.
# Based on Manual X (United Nations, 1986). A good explanation in https://demographicestimation.iussp.org/content/indirect-estimation-child-mortality
# Possibility to use Trussell and Palloni-Helligman variants too
# author: IW
im_brass <- function (CEB, CD, W,
                      age = seq(15, 45, 5),
                      variant = "trussell",
                      mlt_model = NULL,
                      Sex = NULL,
                      date_obs = 0,
                      mam = NULL){
  # panama example in manual x
  # panama <- data.frame(
  #   age = seq(15, 45, 5),
  #   sex = rep("male", 7),
  #   W = c(2695, 2095, 1828, 1605, 1362, 1128, 930),
  #   CEB = c(278, 1380, 2395, 3097, 3444, 3274, 2682),
  #   CD = c(24, 77, 172, 236, 348, 394, 354))
  # W = panama$W
  # CEB = panama$CEB
  # CD = panama$CD
  # age = panama$age
  # Sex = "male"
  # mam = 27

  if(is.null(Sex)) {Sex <- "male"; message("male assumed")}
  if(is.null(mlt_model)) {mlt_model <- "West"; message("west assumed")}
  coef_trussell <- data.frame(
    model = c(rep("North", 7), rep("South", 7), rep("East", 7), rep("West", 7)),
    age = rep(c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"), 4),
    ai = c(1.1119, 1.239, 1.1884, 1.2046, 1.2586, 1.224, 1.1772, 1.0819, 1.2846, 1.2223, 1.1905, 1.1911, 1.1564, 1.1307, 1.1461, 1.2231, 1.1593, 1.1404, 1.154, 1.1336, 1.1201,1.1415, 1.2563, 1.1851, 1.172, 1.1865, 1.1746, 1.1639),
    bi = c(-2.9287,-0.6865, 0.0421, 0.3037, 0.4236, 0.4222, 0.3486, -3.0005, -0.6181, 0.0851,0.2631, 0.3152, 0.3017, 0.2596, -2.2536, -0.4301, 0.0581, 0.1991, 0.2511, 0.2556, 0.2362, -2.707, -0.5381, 0.0633, 0.2341, 0.308, 0.3314, 0.319),
    ci = c(0.8507, -0.2745, -0.5156, -0.5656, -0.5898, -0.5456, -0.4624, 0.8689, -0.3024, -0.4704, -0.4487, -0.4291, -0.3958, -0.3538, 0.6259, -0.2245, -0.3479, -0.3487, -0.3506, -0.3428, -0.3268, 0.7663, -0.2637, -0.4177, -0.4272, -0.4452, -0.4537, -0.4435),
    di = c(1.0921, 1.3207, 1.5996, 2.0779, 2.7705, 4.152, 6.965, 1.09, 1.3079, 1.5173, 1.9399, 2.6157, 4.0794, 7.1796, 1.0959, 1.2921, 1.5021, 1.9347, 2.6197, 4.1317, 7.3657, 1.097, 1.3062, 1.5305, 1.9991, 2.7632, 4.3468, 7.5242),
    ei = c(5.4732, 5.3751, 2.6268, -1.7908, -7.3403, -12.2448, -13.916, 5.4443, 5.5568, 2.6755, -2.2739, -8.4819, -13.8308, -15.388, 5.5864, 5.5897, 2.4692, -2.6419, -8.9693, -14.355, -15.8083, 5.5628, 5.5677, 2.5528, -2.4261, -8.4065, -13.2436, -14.2013),
    fi = c(-1.9672, 0.2133, 4.3701, 9.4126, 14.9352, 19.2349, 19.9542, -1.9721, 0.2021, 4.7471, 10.3876, 16.5153, 21.1866, 21.7892, -1.9949, 0.3631, 5.0927, 10.8533, 17.0981, 21.8247, 22.3005, -1.9956, 0.2962,4.8962, 10.4282, 16.1787, 20.199, 20.0162)
  )
  coef_palloni_helligman <- data.frame(model = c(rep("Latin", 7), rep("Chilean", 7), rep("South_Asian", 7), rep("Far_East_Asian", 7), rep("General", 7)),
                                       age = rep(c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49"), 5),
                                       ai = c(0.6892, 1.3625, 1.0877, 0.75, 0.5605, 0.5024, 0.5326, 0.8274, 1.3129, 1.0632, 0.8236, 0.6895, 0.6098, 0.5615, 0.6749, 1.3716, 1.0899, 0.7694, 0.6156, 0.6077, 0.6952, 0.7194, 1.2671, 1.0668, 0.7833, 0.5765, 0.4115, 0.3071, 0.721, 1.3115, 1.0768, 0.7682, 0.5769, 0.4845, 0.476),
                                       bi = c(-1.6937, -0.3778, 0.0197, 0.0532, 0.0222, 0.0028, 0.0052, -1.5854, -0.2457, 0.0196, 0.0293, 0.0068, -0.0014, 0.004, -1.758, -0.3652, 0.0299, 0.0548, 0.0231, 0.004, 0.0018, -1.3143, -0.2996, 0.0017, 0.0307, 0.0068, 0.0014, 0.0111, -1.4686, -0.336, 0.0109, 0.0439, 0.0176, 0.0034, 0.0071),
                                       ci = c(0.6464, -0.2892, -0.2986, -0.1106, 0.017, 0.0048, 0.0256, 0.5949, -0.2329, -0.1996, -0.0684, 0.0032, 0.0166, 0.0073, 0.6805, -0.2966, -0.2887, -0.0934, 0.0298, 0.0573, 0.0306, 0.5432, -0.2105, -0.2424, -0.1103, -0.0202, 0.0083, 0.0129, 0.5746, -0.2475, -0.2695, -0.109, 0.0038, 0.0036, 0.0246),
                                       di = c(0.0106, -0.0041, 0.0024, 0.0115, 0.0171, 0.018, 0.0168, 0.0097, -0.0031, 0.0021, 0.0081, 0.0119, 0.0141, 0.0159, 0.0109, -0.0041, 0.0024, 0.0108, 0.0149, 0.0141, 0.0109, 0.0093, -0.0029, 0.0019, 0.0098, 0.0165, 0.0213, 0.0251, 0.0095, -0.0034, 0.0021, 0.0105, 0.0165, 0.0187, 0.0189),
                                       ei = c(1.1703, 1.6955, 1.8296, 2.1783, 2.8836, 4.458, 6.9351, 1.3092, 1.6897, 1.8368, 2.2036, 2.9955, 4.7734, 7.4495, 1.1922, 1.7173, 1.8631, 2.1808, 2.7654, 4.1378, 6.4885, 1.2779, 1.7471, 1.9107, 2.3172, 3.2087, 5.1141, 7.6383, 1.2136, 1.7025, 1.836, 2.1882, 2.9682, 4.6526, 7.1425),
                                       fi = c(0.5129, 4.132, 2.902, -2.5688, -10.3282, -17.1809, -19.3871, 1.9474, 4.6176, 2.637, -3.352, -11.4013, -17.885, -19.0513, 0.794, 4.3117, 2.8767, -2.7219, -10.8808, -18.6219, -22.2001, 1.5714, 4.2638, 2.7285, -2.6259, -9.8891, -15.3263, -15.5739, 0.974, 4.1569, 2.8632, -2.6521, -10.3053, -16.692, -18.3021),
                                       gi = c(-0.385, -0.1635, 3.4707, 9.0883, 15.4301, 20.4296, 23.4007, -0.7982, -0.0173, 4.0305, 9.9233, 16.3441, 20.8883, 23.0529, -0.5425, -0.1653, 3.5848, 9.3705, 16.2255, 22.239, 26.4911, -0.6994, -0.0752, 3.5881, 9.0238, 14.7339, 18.2507, 19.7669, -0.5247, -0.1232, 3.522, 9.1961, 15.3161, 19.8534, 22.4168)
  )

  # parity
  Px <- CEB/W
  # death proportion
  Dx <- CD/CEB
  # multipliers
  if(variant == "trussell"){
    if(!mlt_model %in% c("North", "South", "West", "East")) {mlt_model = "West"; message("model was assumed West")}
    k <- as.matrix(coef_trussell[coef_trussell$model == mlt_model, c("ai", "bi", "ci")]) %*% c(1, Px[1]/Px[2], Px[2]/Px[3])
    t <- as.matrix(coef_trussell[coef_trussell$model == mlt_model, c("di", "ei", "fi")]) %*% c(1, Px[1]/Px[2], Px[2]/Px[3])
  }else if(variant == "palloni_helligman"){
    if(is.null(mam)) {mam = 27; message("mam was assumed at 27")}
    if(mlt_model %in% c("North", "South", "West", "East")) {mlt_model = "General"; message("model was assumed General")}
    k <- as.matrix(coef_palloni_helligman[coef_palloni_helligman$model == mlt_model, c("ai", "bi", "ci", "di")]) %*% c(1, Px[1]/Px[2], Px[2]/Px[3], mam)
    t <- as.matrix(coef_palloni_helligman[coef_palloni_helligman$model == mlt_model, c("ei", "fi", "gi")]) %*% c(1, Px[1]/Px[2], Px[2]/Px[3])
  }else{
    stop("variant should be trussell or palloni_helligman")
  }

  # q estimate
  qx <- k * Dx
  lx <- data.frame(age = age, lx = 1 - qx, x = c(1, 2, 3, 5, 10, 15, 20), t, date_t = date_obs - t)
  # what mlt to use
  mlt <- subset(DemoToolsData::modelLTx1, family == mlt_model & sex == Sex) %>%
    dplyr::mutate(lx = lx1/lx1[1], .by = c(source, family, sex, e0))
  # interpolate for each census age
  mlts_interp <- purrr::map_df(1:nrow(lx), function(i){
    la <- lx[i,]
    distances <- mlt$lx[mlt$age==la$x]-la$lx
    best_dist <- distances[abs(distances) %in% sort(abs(mlt$lx[mlt$age==la$x]-la$lx))[1:2]]
    w <- abs(best_dist)/sum(abs(best_dist))
    levels <- (mlt[mlt$age == la$x, ] %>% mutate(dist = lx-la$lx) %>% filter(dist %in% best_dist))$e0
    mlt_interp <- mlt %>%
      summarise(mx_mlt = mx1[e0==levels[1]]*(1-w[1]) + mx1[e0==levels[2]]*w[1],
                qx_mlt = qx1[e0==levels[1]]*(1-w[1]) + qx1[e0==levels[2]]*w[1],
                lx_mlt = lx[e0==levels[1]]*(1-w[1])  + lx[e0==levels[2]]*w[1],
                ex_mlt = ex1[e0==levels[1]]*(1-w[1]) + ex1[e0==levels[2]]*w[1], .by = age) %>%
      rename(age_mlt = age) %>%
      mutate(sex = Sex, x = la$x, lx = la$lx, date_t = la$date_t, .before = 1) %>%
      mutate(age = la$age, .before = 1)
    return(mlt_interp)
  })
  return(list(mi_brass = lx, mlts_interp = mlts_interp))
}
