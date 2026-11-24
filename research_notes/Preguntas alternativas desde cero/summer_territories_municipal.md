# ¿Dónde duermen los residentes durante las vacaciones?

## Dictamen

**Es viable a escala nacional para 2.569 zonas MITMA, con identificación municipal inequívoca para 1.643 de los 8.131 códigos INE del cruce.** Estos municipios concentran el 88,8 % de la población de referencia del fichero de zonas. Hay 924 zonas `_AM` y un caso adicional de lindes incompatibles con una unión directa: Daya Vieja (03903) aparece tanto en su zona homónima como en Almoradí (03015). La precisión municipal se pierde también en el otro extremo de cada enlace si el destino es una agrupación.

El conjunto incluye julio y agosto completos, el verano combinado, Semana Santa, Todos los Santos y una referencia laborable de febrero. Agosto de 2024 tiene masa positiva en las 2.569 zonas de residencia y 2.568 zonas de destino. `19103_AM` carece de masa publicada como destino en agosto. Todos los municipios con correspondencia inequívoca tienen estimación de destino en agosto. El 92,7 % de las combinaciones teóricas zona de origen × zona de destino × día de agosto carecen de fila. Esa cifra incluye desplazamientos inexistentes y ausencias por causas no identificables en la tabla. Los enlaces con poca masa y las entropías merecen más cautela que los agregados de presencia. La coincidencia con el INE no permite reconstruir celdas ausentes. [Metodología MITMA](https://cdn.transportes.gob.es/portal-web-transportes/ministerio/bigdata/a3_informe_metodologico_resumen_v3.pdf).

## Estabilidad y contraste

- Del 12 al 18 de agosto, el Spearman del multiplicador zonal frente a 2024 es aproximadamente 0,98 en 2022 y 2023. El primer destino coincide en el 72 % y 76 % de las zonas comparables, respectivamente. La comparación mantiene fijo el denominador territorial y mezcla calendarios semanales diferentes. La estabilidad de los multiplicadores no implica estabilidad de todos los enlaces ni de las colonias.
- En agosto de 2024, 84.794 de los 85.952 pares positivos del INE dentro del universo de municipios directos tienen también masa MITMA positiva. En los 1.441 orígenes con cinco o más pares compartidos, la mediana de Spearman de los recuentos es 0,577 y la coincidencia mediana es 3 de 5 destinos, restringidos a los pares comunes. El INE aporta una observación mensual de turismo interprovincial con reglas y unidades diferentes. Puede sostener un producto municipal alternativo para los pares publicados, con cobertura limitada por viaje, muestra y difusión. No permite estimar un fin de semana concreto ni incluye a los residentes que duermen en casa para reconstruir el multiplicador. [Libro INE de 2024](https://www.ine.es/experimental/turismo_moviles/exp_tmov_interno_mun_2024.xlsx) · [Metodología INE](https://ine.es/experimental/turismo_moviles/exp_moviles_turismo_interno.pdf).

## Valor para Arraigo

En la ficha territorial, mostrar «Presencia nocturna estacional» permitiría distinguir un lugar con actividad residencial estable de otro que cambia sensiblemente en agosto o en festivos. «Dónde duermen los residentes fuera de su zona» aportaría contexto sobre conexiones de fin de semana y verano. Una proporción por origen y el destino dominante son más fáciles de interpretar que una razón de gravedad. La presencia aporta contexto para valorar la variación estacional de la actividad local. La presión en servicios, el alquiler y el motivo de los viajes requieren sus propias fuentes.

## Contrato propuesto

La inspección en solo lectura confirmó que Arraigo distribuye JSON comprimido bajo `front/public/data/v1/` y documenta fuente, unidad y disponibilidad por escala en sus catálogos. Propuesta de artefactos futuros:

| Ruta respecto a `front/public/data/v1/` | Campos |
|---|---|
| `context/overnight_seasonality_municipal_2024.json.gz` | `ine_code`, `zone_id`, `period`, `days`, `origin_reporting_days`, `destination_reporting_days`, `present_night_mean`, `resident_population`, `population_multiplier`, `present_vs_ordinary_ratio`, `inbound_night_mean`, `dominant_origin_zone_id`, `dominant_origin_share`, `origin_entropy`, `away_share`, `coastal_share`, `own_province_share`, `mean_distance_km`, `destination_entropy`. Clave `ine_code,period`, con código INE de cinco dígitos en texto. |
| `context/overnight_top_destinations_2024.json.gz` | `origin_ine_code`, `period`, `rank`, `destination_zone_id`, `destination_ine_code` (nullable), `destination_name`, `destination_municipality_count`, `night_mean`, `away_share`, `distance_km`. |
| `context/overnight_mitma_coverage_2024.json.gz` | `ine_code`, `zone_id`, `granularity` (`municipality`, `aggregate_zone`, `split_municipality`), `zone_name`, `municipality_count`. Conserva una fila por correspondencia, incluida la múltiple. |
| `metadata/sources_catalog.json` | Fuente `mitma_overnight_stays_v2`, registrando `source_id,name,publisher,source_url,format,data_date,spatial_levels,native_granularity,variables,terms_notes`. |
| `metadata/variables_metadata.json` | Entradas `overnight_population_multiplier` y `overnight_present_vs_ordinary_ratio`, con `indicator_id,name,category_id,unit,source_ids,aggregation_method,notes,availability_by_level`. Disponibilidad municipal selectiva y `scoring_enabled: false`. |

El cociente frente a una noche ordinaria sería `present_night_mean(period)/present_night_mean(ordinary_2024)`, disponible si la cobertura diaria es completa en ambos periodos. Para el multiplicador, se recomienda comprobar el año del denominador poblacional antes de mostrarlo como padrón. La ficha de un municipio agrupado puede explicar la cobertura y enlazar a la agrupación con su nombre, sin atribuirle el valor del grupo. Los datos de colonia quedan en el conjunto de investigación.
