# TFG
En este repositorio se muestran los datos empleados, los códigos desarrollados, tanto con emcee como con Cobaya, así como los códigos empleados para obtener los coner plots como las tablas, en mi TFG sobre "Inferencia bayesiana en modelos de energía oscura con datos cosmológicos"

## Directorios:
- Modelos: distinguimos entre emcee y Cobaya
- emcee: encontramos los 4 modelos desarrollados, además de un repositorio con los códigos para obtener los corner plots de dichos modelos.
- Cobaya: encontramos los 2 modelos desarrollados, además de un repositorio con los códigos para obtener las tablas de dichos modelos, ya que la interfaz de getdist permite realizar los corner plots con los archivos de Cobaya fácilmente.

## Modelos analizados
- **ΛCDM flat**: [h, Ωₘ, M]
- **ΛCDM**: [h, Ωₘ, Ω_k, M]
- **wCDM**: [h, Ωₘ, M, w]
- **w₀wₐCDM**: [h, Ωₘ, M, w₀, wₐ]

## Metodología
- Datos de Pantheon+SH0ES para supernovas Ia
- MCMC samples a partir de [`emcee`] (https://emcee.readthedocs.io) y de [`Cobaya`] (https://cobaya.readthedocs.io/en/latest/index.html)
- Likelihoods basada en la magnitud aparente de supernovas Ia
- Corner plots a partir de de getdist para las distribuciones de los parámetros

## Requisitos
- Python => 3.9
- numpy, matplotlib, scipy, emcee, getdist, cobaya

## Autor
TFG -- Física
Universidad de Córdoba
Ismael Antequera García

## License
MIT
