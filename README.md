# Simulador del Dogma Central de la Biología Molecular

**Bioinformática — Práctica 1**

## Autores

- Gabriel Munteanu Cabrera
- Andrea Dumpiérrez Medina

**ULPGC – GCID (Ciencias e Ingeniería de Datos)**

## Descripción

Este proyecto consiste en el desarrollo de un simulador en Python que representa los principales procesos del dogma central de la biología molecular: replicación del ADN, transcripción a ARN mensajero y traducción para obtener una secuencia de aminoácidos.

El programa permite seguir cada etapa mediante una interfaz de texto, mostrando las moléculas, enzimas y transformaciones que intervienen en los procesos.

## Funcionalidades

- **Entrada del ADN:** permite introducir una secuencia manualmente, cargarla desde un archivo FASTA o generarla aleatoriamente.
- **Replicación:** representa la separación de las hebras de ADN, la formación de la cadena líder y la cadena rezagada mediante fragmentos de Okazaki, así como la participación de cebadores y enzimas.
- **Transcripción:** genera una secuencia de ARN mensajero a partir de la cadena molde de ADN, aplicando las reglas de complementariedad.
- **Traducción:** identifica el codón de inicio AUG, traduce los codones a aminoácidos y finaliza al encontrar un codón STOP.
- **Visualización:** muestra los resultados de cada etapa y un resumen final del flujo de información genética.

## Archivos del proyecto

- `simulador.py`: código principal del simulador.
- `lacZ_Ecoli_NC_000913.3.fasta`: secuencia del gen *lacZ* de *Escherichia coli*, utilizada para probar el programa.

## Requisitos

- Python 3.

No es necesario instalar bibliotecas externas.

## Ejecución

Desde una terminal, situarse en la carpeta del proyecto y ejecutar:

```bash
python simulador.py
```
El programa mostrará tres opciones:

1. Introducir ADN manualmente.
2. Leer ADN desde un archivo FASTA.
3. Generar ADN aleatoriamente.

Si se selecciona la segunda opción, puede utilizarse el archivo `lacZ_Ecoli_NC_000913.3.fasta` incluido en el repositorio.

## Consideraciones

Las secuencias largas se muestran de forma abreviada para mejorar la legibilidad de los resultados.

Algunos aspectos, como el tamaño de los cebadores y de los fragmentos de Okazaki, se han simplificado para facilitar la representación de los procesos.
