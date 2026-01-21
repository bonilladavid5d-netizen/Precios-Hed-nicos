###############################################
# Trabajo final 
###############################################


library(sf)
library(dplyr)

setwd("/Users/davidbonilla/Downloads/Final_espacial")

manzanas  <- st_read("MGN_ANM_MANZANA.shp")
valor_ref <- st_read("Valor_Ref_M_2024.shp")

if (st_crs(valor_ref) != st_crs(manzanas)) {
  valor_ref <- st_transform(valor_ref, st_crs(manzanas))
}

sf_use_s2(FALSE)

manzanas_bog <- manzanas[valor_ref, ]
manzanas_bog_valor <- st_join(manzanas_bog, valor_ref, join = st_intersects)

manzanas_bog_valor <- manzanas_bog_valor |>
  mutate(ln_valor_m2 = log(ValRef))

est_tm <- st_read("Estación_troncal.shp")

crs_metros <- 3116

manzanas_bog_valor_m <- st_transform(manzanas_bog_valor, crs_metros)
est_tm_m              <- st_transform(est_tm,           crs_metros)

idx_tm  <- st_nearest_feature(manzanas_bog_valor_m, est_tm_m)

dist_tm <- st_distance(manzanas_bog_valor_m, est_tm_m[idx_tm, ], by_element = TRUE)

manzanas_bog_valor_m <- manzanas_bog_valor_m |>
  mutate(dist_tm = as.numeric(dist_tm))

aeropuerto_wgs <- st_sfc(
  st_point(c(-74.1469, 4.7016)),
  crs = 4326
)

aeropuerto_m <- st_transform(aeropuerto_wgs, st_crs(manzanas_bog_valor_m))

manz_centroid_m <- st_centroid(manzanas_bog_valor_m)

dist_mat  <- st_distance(manz_centroid_m, aeropuerto_m)
dist_aero <- as.numeric(dist_mat)

manzanas_bog_valor_m <- manzanas_bog_valor_m |>
  mutate(
    dist_aero    = dist_aero,
    dist_aero_km = dist_aero / 1000
  )

centro_wgs <- st_sfc(
  st_point(c(-74.08175, 4.60971)),
  crs = 4326
)

centro_m <- st_transform(centro_wgs, st_crs(manzanas_bog_valor_m))

dist_centro_m  <- st_distance(manz_centroid_m, centro_m)
dist_centro    <- as.numeric(dist_centro_m)

manzanas_bog_valor_m <- manzanas_bog_valor_m |>
  mutate(
    dist_centro    = dist_centro,
    dist_centro_km = dist_centro / 1000
  )

manzanas_bog_valor_m <- manzanas_bog_valor_m |>
  mutate(
    log_densidad = log1p(DENSIDAD),
    share_secund   = ifelse(TP27_PERSO > 0, TP51SECUND / TP27_PERSO, NA_real_),
    share_superior = ifelse(TP27_PERSO > 0, TP51SUPERI / TP27_PERSO, NA_real_)
  )

predio <- st_read("bd_uaecd.gpkg", layer = "Predio")

controles_loc <- predio |>
  st_drop_geometry() |>
  group_by(PreCodLoc) |>
  summarise(
    estrato_loc_prom      = mean(PreEstrato[PreEstrato > 0], na.rm = TRUE),
    area_const_loc_prom   = mean(PreAConst,           na.rm = TRUE),
    antiguedad_loc_prom   = mean(2025 - as.numeric(PreVForma), na.rm = TRUE),
    .groups = "drop"
  )

manzanas_bog_valor_m <- manzanas_bog_valor_m |>
  mutate(CD_LC_CM = as.numeric(CD_LC_CM)) |>
  left_join(
    controles_loc,
    by = c("CD_LC_CM" = "PreCodLoc")
  )

base_modelo <- manzanas_bog_valor_m |>
  st_drop_geometry() |>
  mutate(
    dist_tm_km   = dist_tm   / 1000,
    dist_aero_km = dist_aero / 1000
  ) |>
  select(
    ln_valor_m2,
    dist_aero_km,
    dist_tm_km,
    log_densidad,
    share_secund,
    share_superior,
    antiguedad_loc_prom,
    dist_centro_km,
    estrato_loc_prom,
    area_const_loc_prom
  ) |>
  filter(
    dist_aero_km <= 15,
    !if_any(everything(), is.na)
  )

library(dplyr)
library(broom)
library(ggplot2)
library(sf)

modelo_aero_final <- lm(
  ln_valor_m2 ~ dist_aero_km +
    dist_tm_km +
    log_densidad +
    share_secund +
    share_superior +
    antiguedad_loc_prom +
    dist_centro_km +
    estrato_loc_prom +
    area_const_loc_prom,
  data = base_modelo
)

summary(modelo_aero_final)

base_modelo <- base_modelo |>
  mutate(
    dist_aero_km_c  = dist_aero_km - mean(dist_aero_km, na.rm = TRUE),
    dist_aero_km_c2 = dist_aero_km_c^2
  )

modelo_aero_quad <- lm(
  ln_valor_m2 ~ dist_aero_km_c + dist_aero_km_c2 +
    dist_tm_km +
    log_densidad +
    share_secund +
    share_superior +
    antiguedad_loc_prom +
    dist_centro_km +
    estrato_loc_prom +
    area_const_loc_prom,
  data = base_modelo
)

summary(modelo_aero_quad)

base_modelo <- base_modelo |>
  mutate(
    log_dist_aero = log(dist_aero_km + 0.1)
  )

modelo_aero_log <- lm(
  ln_valor_m2 ~ log_dist_aero +
    dist_tm_km +
    log_densidad +
    share_secund +
    share_superior +
    antiguedad_loc_prom +
    dist_centro_km +
    estrato_loc_prom +
    area_const_loc_prom,
  data = base_modelo
)

summary(modelo_aero_log)

base_modelo <- base_modelo |>
  mutate(
    anillo_aero = cut(
      dist_aero_km,
      breaks = c(-Inf, 2, 5, 10, Inf),
      labels = c("0-2km", "2-5km", "5-10km", ">10km")
    ),
    anillo_aero = stats::relevel(anillo_aero, ref = ">10km")
  )

modelo_aero_anillos <- lm(
  ln_valor_m2 ~ anillo_aero +
    dist_tm_km +
    log_densidad +
    share_secund +
    share_superior +
    antiguedad_loc_prom +
    dist_centro_km +
    estrato_loc_prom +
    area_const_loc_prom,
  data = base_modelo
)

summary(modelo_aero_anillos)

coefs <- tidy(modelo_band, conf.int = TRUE)

coefs_dist <- subset(coefs, grepl("^dist_aero_band", term))

coefs_dist$band <- gsub("dist_aero_band", "", coefs_dist$term)

coefs_dist$band <- factor(
  coefs_dist$band,
  levels = c("0-2", "2-4", "4-6", "6-8", "8-15")
)

ggplot(coefs_dist, aes(x = band, y = estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Banda de distancia al aeropuerto (km)",
    y = "Efecto sobre ln(precio m²) vs. 8–15 km",
    title = "Efecto de vivir cerca del aeropuerto sobre el precio por m²"
  )

manzanas_10km <- manzanas_bog_valor_m |>
  dplyr::filter(dist_aero_km <= 10)

ggplot() +
  geom_sf(data = manzanas_10km, aes(fill = ValRef), color = NA) +
  labs(
    title = "Valor m² en un radio de 10 km alrededor del aeropuerto",
    fill  = "valor m²"
  ) +
  theme_minimal()

install.packages(c("broom", "flextable", "officer"))

library(broom)
library(dplyr)
library(flextable)
library(officer)


export_lm_to_word <- function(model, titulo = "Modelo", file = "modelo.docx") {
  
  tab <- broom::tidy(model) |>
    mutate(
      # Asteriscos de significancia
      Signif = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.10  ~ ".",
        TRUE            ~ ""
      ),
      estimate   = round(estimate,   4),
      std.error  = round(std.error,  4),
      statistic  = round(statistic,  3),
      p.value    = round(p.value,    4)
    ) |>
    dplyr::select(
      Variable     = term,
      Coeficiente  = estimate,
      `Error estándar` = std.error,
      t            = statistic,
      `p-valor`    = p.value,
      Signif
    )
  
  ft <- flextable(tab)
  ft <- autofit(ft)
  ft <- add_header_row(ft, values = titulo, colwidths = ncol(tab))
  
  doc <- read_docx()
  doc <- body_add_flextable(doc, ft)
  print(doc, target = file)
}
export_lm_to_word(
  modelo_aero_final,
  titulo = "Modelo lineal distancia",
  file   = "modelo_aero_final.docx"
)

export_lm_to_word(
  modelo_aero_quad,
  titulo = "Modelo cuadrático distancia",
  file   = "modelo_aero_quad.docx"
)

export_lm_to_word(
  modelo_aero_log,
  titulo = "Modelo log(distancia)",
  file   = "modelo_aero_log.docx"
)

export_lm_to_word(
  modelo_aero_anillos,
  titulo = "Modelo por anillos de distancia",
  file   = "modelo_aero_anillos.docx"
)

