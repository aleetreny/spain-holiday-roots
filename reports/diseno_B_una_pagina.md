# Pregunta B: diseño "primero la red" (1 página)

## Pregunta
- **Título + pregunta (EN, elegido por el usuario):** ***Roots or beaches?*** *On holiday, do Spaniards go back to the villages their families left?*
- **ES:** ¿Vuelve España en vacaciones a los pueblos que dejaron sus familias?
- **Contexto en una línea (primer párrafo del informe):** *Between 1950 and 1975 millions of Spaniards left their villages for Madrid, Barcelona, the Basque Country and other cities — the "rural exodus" — and many villages lost most of their people.* [cifra exacta y cita: por verificar]
- **Versión operativa (sección 2):** comparado con una noche normal, ¿la red de pernoctaciones en vacaciones se reorganiza hacia los lugares donde nacieron las familias de los residentes, más de lo que explican la distancia, la población, la renta y la playa?
- **Respuesta esperada:** sí, sobre todo en verano y en el interior; no es solo turismo (quitar la costa refuerza el efecto). Matiz honesto: la migración reciente también pesa, sobre todo en Navidad.
- **Por qué no es obvia:** España es el 2.º destino turístico del mundo; lo esperable es que playa, renta y segunda residencia dominen el mapa de vacaciones.

## Redes (mismos 2.484 nodos = zonas MITMA de la España peninsular)
| Red | Nodo | Enlace | Papel |
|---|---|---|---|
| Noche normal (miér. 14-feb-2024) | zona | residentes de i que duermen en j (≥10) | referencia, densidad 0,86 % |
| **Noche de vacaciones (15-ago-2024)** | zona | igual | red principal, densidad 3,3 % |
| Todos los Santos (1-nov-2024) | zona | igual | control sin playa, 1,6 % |
| Raíces | provincia / 137 ciudades | nacidos en p que viven en i (INE) | red de referencia para comparar (no se le piden propiedades básicas: es densa) |
Navidad, Semana Santa, media de agosto y años 2022–23 → apéndice.

**Representatividad (decidido con el usuario: las 3 noches valen si son robustas).** Cada noche principal se repite con sus noches "hermanas" y se da la mediana y el rango, no un único valor:
- 15-ago → cada noche del 10 al 18 de agosto, por separado.
- 14-feb → todos los martes, miércoles y jueves de febrero.
- 1-nov → 31-oct y 2-nov.
- Las tres → 2022 y 2023 (el 1-nov-2023 falta por un corte del MITMA; se usa 2022).
- Criterio: el resultado vale si todas las noches hermanas dan el mismo signo y la noche principal no es la extrema del grupo.

## Métricas por apartado de la memoria
- **Metodología (propiedades básicas, cada red vs su versión aleatorizada):** N, E, densidad, componente gigante; distribución de grado (CCDF de k_in y k_out); camino medio y diámetro; clustering vs red recableada; reciprocidad.
- **Resultado 1, la red cambia de forma:** asortatividad r (real vs recableada → está forzada por los grados, Zhou–Mondragón); club de ricos normalizado ρ(k) (1,25 → 0,99); k-core de los pueblos vaciados (comparando pueblos del mismo tamaño); BOSAM ordenada por población.
- **Resultado 2, la red de vacaciones se parece a la de raíces:** QAP con modelo nulo de gravedad espacial, QAP parcial controlando la noche normal, diseño de ciudades vecinas (misma provincia, distintas raíces); PPML con control de renta como confirmación.
- **Resultado 3, territorios de verano:** comunidades bipartitas (ciudades × destinos) con modularidad de Barber bajo BiCM y bajo gravedad espacial; cuántas comunidades son "archipiélagos" no contiguos; enriquecimiento de raíces dentro de cada territorio.

## Modelos nulos (solo tres familias)
1. **Recableado que conserva el grado** (25×): asortatividad, club de ricos, k-core, clustering.
2. **Gravedad espacial doblemente restringida**: quita distancia y tamaño (QAP y comunidades).
3. **Permutaciones QAP** (2.000; también solo entre ciudades de la misma provincia) y **BiCM** para la red bipartita.
Barber ponderado, JDD y Monte Carlo de contigüidad → apéndice.

## Figuras (5) y tablas (2)
1. **Mapa de las tres redes** (noche normal, 15 de agosto, raíces) con los mismos nodos. *Falta.*
2. **Cambio de forma:** CCDF de grado, r real vs recableada por fecha, curvas ρ(k), k-core por quintil de despoblación. (rehacer `net_B_topology_shift.png`)
3. **BOSAM** febrero vs agosto (la "cortina"). (ya existe como panel; separar)
4. **Raíces:** dispersión ciudades vecinas (ej. Barcelona → Andalucía) + distribución nula del QAP. (rehacer a partir de `net_B_qap_and_territories.png` y `fig9`)
5. **Mapa de territorios de verano**: BiCM vs nulo espacial.
- Tabla 1: propiedades básicas de cada red y de su versión aleatorizada. Tabla 2: resumen de los tres resultados con su nulo, z y p.
