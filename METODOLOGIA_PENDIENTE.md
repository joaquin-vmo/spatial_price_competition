# Pendientes metodológicos

Estado al 31 de julio de 2026. Lista de lo que falta, ordenada por el valor que
agregaría a la defensa del trabajo frente al costo de hacerlo.

---

## 1. Sin solución con los datos actuales

### 1.1 Placebo sobre demanda local

**Qué es.** El contraste que verifica que la entrada no coincide con un aumento
de la demanda local que también movería los precios. Fischer, Martin y
Schmidt-Dengler lo hacen con dos fuentes: 2.000 contadores de flujo vehicular y
una grilla poblacional de 1×1 km. Muestran que ni el tráfico ni la población se
mueven alrededor de la entrada.

**Por qué no se puede.** Las ventas disponibles son nacionales-mensuales
(`monthly_sales_1.xlsx`) y regionales-anuales (`monthly_sales_2.xlsx`). Ninguna
identifica demanda al nivel del mercado local de 2,5 km.

**Sustituto posible y su límite.** Parque automotor y proyecciones de población
comunal del INE. La comuna sigue siendo una unidad enorme frente a un radio de
2,5 km: comunas como Santiago o Puente Alto contienen varios mercados locales
distintos. La recomendación es declararlo ausente en limitaciones antes que
aproximarlo mal.

**Estado actual.** El hazard de timing (`12_atenuacion_y_timing.R`) cubre
parcialmente el flanco: muestra que lo que predice la entrada es la densidad de
estaciones y no la rentabilidad local. Pero no puede descartar que la
densificación venga acompañada de crecimiento de demanda. Está declarado como
limitación en la subsección de validez.

### 1.2 Fechas de permiso (DOM, SEC)

**Qué es.** Redefinir el tratamiento en la fecha de inicio de construcción o de
otorgamiento del permiso, en vez de la fecha de apertura. Convertiría el
argumento conceptual del rezago administrativo en evidencia.

**Por qué no se puede.** No existe registro público consolidado. Habría que
pedirlo por transparencia a cada Dirección de Obras Municipales o a la SEC.

**Mitigación en uso.** El test de anticipación semanal
(`20_anticipacion_seleccion.R`) muestra trayectorias planas en las 25 semanas
previas a la apertura, lo que hace poco probable que el efecto haya empezado
antes.

---

## 2. Factible, no hecho

### 2.1 Efectos distribucionales por cuantiles (RIF / Chernozhukov et al.)

**Prioridad: alta.** Es lo que separa a Fischer de un event study convencional y
es el aporte que más levantaría la tesis.

El script `16_mercado_max_media_min.R` ya entrega la intuición: el mínimo del
mercado cae 5,5 $/L y el máximo no se mueve, a 1 km. Falta la versión formal:

- Regresiones RIF (`firpo.pdf` ya está en `literature/`) para efectos por
  cuantil de la distribución de precios.
- Test de dominancia estocástica de primer orden (Tabla 3 de Fischer).
- El *value of information*: diferencia entre precio medio y mínimo del mercado
  como medida de la ganancia de buscar.

**Costo:** un script. Los insumos están todos construidos.

### 2.2 Callaway y Sant'Anna

**Prioridad: media.** El paper está en `literature/`. Sun-Abraham ya está
implementado (`19_staggered_persistencia.R`) y coincide con TWFE, así que
agregar un segundo estimador robusto es marginal. Vale como robustez adicional
si un lector lo pide.

### 2.3 Mapa de dispersión geográfica y temporal

**Prioridad: media-baja, costo bajo.** La Figura 1 de Fischer: las 403 entradas
focales sobre el mapa de Chile más el histograma temporal. Sustenta el argumento
de que un shock confusor tendría que replicar esa dispersión para invalidar el
diseño.

### 2.4 Balance ampliado con amenidades

**Prioridad: baja, costo muy bajo.** El vector de servicios de la CNE —tienda,
lavado, baño, farmacia, surtidor de camiones— está en la base y no se usa en la
tabla de balance. Advertencia: solo se dispone de la última actualización de
cada estación, no de su evolución, así que sirve para balance pero no como
control variable en el tiempo.

### 2.5 Diseño de pass-through como supuesto más débil

**Prioridad: media, costo alto.** Identificar un *cambio en la respuesta* a los
shocks MEPCO alrededor de la entrada requiere únicamente que el timing no se
correlacione con cambios en la elasticidad de respuesta, no con el nivel de
demanda local. Es un supuesto estrictamente más débil que el actual.

Los scripts `05_local_projection.R`, `06_distributed_lags.R` y
`07_passthrough_analysis.R` ya estiman pass-through por número de competidores.
Falta reindexarlos a tiempo-evento respecto de la entrada.

### 2.6 Heterogeneidad por marca de la entrante

**Prioridad: baja, costo bajo.** La bandera blanca es el principal entrante (184
de 482) y la FNE sostiene que inyecta mayor presión competitiva. Contrastarlo
sería un resultado con lectura de política directa.

---

## 3. Decisiones abiertas

### 3.1 Radio de la especificación principal

La tesis usa 2,5 km en todo el análisis y adopta de Fischer solo el formato de
reporte (bins de seis meses, ventana ±4). El script `17_formato_fischer.R`
produce además la versión a 1 y 2 km, que son los radios de Fischer.

Si se prefiriera adoptar 1 o 2 km como especificación principal, habría que
rehacer al mismo radio la atenuación por distancia (`12`) y la heterogeneidad
por margen previo (`14`) para que la sección sea internamente coherente.

### 3.2 Gasolina de 97 octanos

Es el combustible más débil de manera consistente: efecto menor, menos preciso, y
en varias especificaciones no significativo. Es el combustible premium, de menor
volumen y con mayor peso de la política de marca. Conviene decidir si se reporta
en el cuerpo con una nota o se relega al apéndice.

---

## 4. Limitaciones de datos ya documentadas

Están escritas en el Anexo B y no requieren acción, pero conviene tenerlas
presentes al responder preguntas.

- **Fechas de salida no confiables.** El cambio de régimen de reporte de la CNE
  concentra 100 bajas en diciembre de 2022. El diseño usa solo entradas, que sí
  son fechables.
- **Censura izquierda.** 1.548 de las 2.030 estaciones ya reportaban en 2012 y su
  fecha de apertura es desconocida.
- **Cohorte 2013 excluida** del tratamiento: 79 entradas, de las cuales 43 son
  bandera blanca, patrón de incorporación tardía al reporte.
- **Duplicados con sufijo `a`.** 31 pares corregidos en el build. Quien retome
  estos datos debe aplicar la misma corrección: sin ella se inflan los conteos
  de competidores y se fabrican salidas.
- **MEPCO es una referencia nacional**, no regional. Fischer dispone de precios
  mayoristas de nueve regiones. Cualquier diferencial de costo regional queda
  dentro del margen medido y solo es absorbido si es constante en el tiempo.
- **Sin datos de cantidades.** No se observan volúmenes por estación, de modo que
  no puede estimarse el efecto sobre ingresos, como sí hace Arcidiacono et al.
