import random


# ============================================================
# 1. TABLAS Y CONSTANTES BIOLÓGICAS
# ============================================================

COMPLEMENTO_ADN = {"A": "T", "T": "A", "C": "G", "G": "C"}

COMPLEMENTO_ARN = {"A": "U", "T": "A", "C": "G", "G": "C"}

BASES_VALIDAS = set("ATCG")

CODON_INICIO = "AUG"

STOP_CODONES = {
    "UAA": "Ocre",
    "UAG": "Ámbar",
    "UGA": "Ópalo"
}

CODIGO_GENETICO = {
    "UUU": "Fenilalanina", "UUC": "Fenilalanina",
    "UUA": "Leucina", "UUG": "Leucina",
    "CUU": "Leucina", "CUC": "Leucina",
    "CUA": "Leucina", "CUG": "Leucina",

    "AUU": "Isoleucina", "AUC": "Isoleucina",
    "AUA": "Isoleucina", "AUG": "Metionina",

    "GUU": "Valina", "GUC": "Valina",
    "GUA": "Valina", "GUG": "Valina",

    "UCU": "Serina", "UCC": "Serina",
    "UCA": "Serina", "UCG": "Serina",

    "CCU": "Prolina", "CCC": "Prolina",
    "CCA": "Prolina", "CCG": "Prolina",

    "ACU": "Treonina", "ACC": "Treonina",
    "ACA": "Treonina", "ACG": "Treonina",

    "GCU": "Alanina", "GCC": "Alanina",
    "GCA": "Alanina", "GCG": "Alanina",

    "UAU": "Tirosina", "UAC": "Tirosina",

    "CAU": "Histidina", "CAC": "Histidina",
    "CAA": "Glutamina", "CAG": "Glutamina",

    "AAU": "Asparagina", "AAC": "Asparagina",
    "AAA": "Lisina", "AAG": "Lisina",

    "GAU": "Ácido aspártico", "GAC": "Ácido aspártico",
    "GAA": "Ácido glutámico", "GAG": "Ácido glutámico",

    "UGU": "Cisteína", "UGC": "Cisteína",
    "UGG": "Triptófano",

    "CGU": "Arginina", "CGC": "Arginina",
    "CGA": "Arginina", "CGG": "Arginina",

    "AGU": "Serina", "AGC": "Serina",
    "AGA": "Arginina", "AGG": "Arginina",

    "GGU": "Glicina", "GGC": "Glicina",
    "GGA": "Glicina", "GGG": "Glicina"
}

# Tamaños reducidos para poder visualizar la simulación.
TAM_CEBADOR = 3
TAM_OKAZAKI = 6
MAX_MOSTRAR = 80 # Número máximo de bases que mostramos completas.

# ============================================================
# 2. FUNCIONES AUXILIARES
# ============================================================

def titulo(texto):
    print("\n" + "=" * 65)
    print(texto)
    print("=" * 65)


def limpiar_adn(secuencia):
    """Elimina espacios y convierte la secuencia a mayúsculas."""
    return "".join(secuencia.split()).upper()


def validar_adn(adn):
    """Comprueba que solo existan A, T, C y G."""
    return len(adn) > 0 and all(base in BASES_VALIDAS for base in adn)


def complementaria(adn):
    """Devuelve la cadena complementaria de ADN."""
    return "".join(COMPLEMENTO_ADN[base] for base in adn)


def complementaria_inversa(adn):
    """Devuelve la complementaria inversa."""
    return complementaria(adn)[::-1]


def abreviar(secuencia):
    """Evita imprimir completas las secuencias muy largas."""
    if len(secuencia) <= MAX_MOSTRAR:
        return secuencia
    mitad = MAX_MOSTRAR // 2
    return secuencia[:mitad] + "..." + secuencia[-mitad:]


def mostrar_adn(codificante, molde, nombre="ADN"):
    """Muestra las dos cadenas antiparalelas."""
    print(f"\n{nombre}:")
    print("5' -", abreviar(codificante), "- 3'")

    if len(codificante) <= MAX_MOSTRAR:
        print("     " + "|" * len(codificante))

    print("3' -", abreviar(molde), "- 5'")
    print("Longitud:", len(codificante), "pb")

# ============================================================
# 3. ENTRADA DEL ADN
# ============================================================

def pedir_adn_manual():
    """Permite introducir una secuencia por teclado."""
    while True:
        adn = input("\nIntroduce la hebra codificante de ADN: ")
        adn = limpiar_adn(adn)

        if validar_adn(adn):
            return adn

        print("ERROR: utiliza únicamente A, T, C y G.")


def leer_fasta(ruta):
    """Lee una secuencia de ADN desde un archivo FASTA."""
    secuencia = []
    cabecera = None

    with open(ruta, "r", encoding="utf-8") as archivo:
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue
            if linea.startswith(">"):
                if cabecera is None:
                    cabecera = linea[1:]
                continue
            secuencia.append(linea)
    adn = limpiar_adn("".join(secuencia))

    if not validar_adn(adn):
        raise ValueError("El FASTA no contiene una secuencia válida de ADN.")
    return adn, cabecera


def generar_adn_aleatorio(num_codones=7):
    """
    Genera ADN con ATG inicial, varios codones intermedios
    y un codón STOP final.
    """
    codones = [
        codon for codon in CODIGO_GENETICO
        if codon != CODON_INICIO ]
    arnm = CODON_INICIO

    for _ in range(num_codones):
        arnm += random.choice(codones)
    arnm += random.choice(list(STOP_CODONES))
    
    return arnm.replace("U", "T")


def elegir_adn():
    """Permite seleccionar el origen del ADN."""
    while True:
        titulo("SELECCIÓN DE LA SECUENCIA")
        print("""
1. Introducir ADN manualmente
2. Leer ADN desde un archivo FASTA
3. Generar ADN aleatoriamente
""")
        opcion = input("Opción [1/2/3]: ").strip()
        if opcion == "1":
            return pedir_adn_manual()
       
        elif opcion == "2":
            ruta = input("Ruta del archivo FASTA: ").strip().strip('"')

            try:
                adn, cabecera = leer_fasta(ruta)
                print("\nFASTA cargado correctamente.")
                if cabecera:
                    print("Cabecera:", cabecera)
                print("Longitud:", len(adn), "nucleótidos")
                return adn
            
            except FileNotFoundError:
                print("\nERROR: no se encuentra el archivo.")

            except (ValueError, OSError) as error:
                print("\nERROR:", error)

        elif opcion == "3":
            adn = generar_adn_aleatorio()
            print("\nADN generado:")
            print(adn)
            return adn

        else:
            print("\nSelecciona 1, 2 o 3.")

# ============================================================
# 4. REPLICACIÓN
# ============================================================

def cebador_arn(fragmento):
    """Genera un pequeño cebador de ARN."""
    longitud = min(TAM_CEBADOR, len(fragmento))
    return fragmento[:longitud].replace("T", "U")


def replicacion(adn):
    """Simula la replicación semiconservativa del ADN."""
    titulo("PASO 1: REPLICACIÓN (ADN -> ADN)")
   
    molde = complementaria(adn)
    mostrar_adn(adn, molde, "ADN ORIGINAL")

    # --------------------------------------------------------
    # Apertura
    # --------------------------------------------------------

    print("\n[HELICASA]")
    print("Separa las dos cadenas rompiendo los puentes de hidrógeno.")

    print("\n[PROTEÍNAS SSB]")
    print("Mantienen separadas y estabilizadas las cadenas de ADN.")

    # --------------------------------------------------------
    # Cadena líder
    # --------------------------------------------------------

    print("\n--- CADENA LÍDER ---")
    print("Se sintetiza continuamente en dirección 5' -> 3'.")

    lider = complementaria(molde)
    cebador = cebador_arn(lider)

    print("\n[PRIMASA]")
    print("Cebador de ARN: 5'-" + cebador + "-3'")

    print("\n[ADN POLIMERASA III]")
    print("Alarga la cadena desde el cebador en dirección 5' -> 3'.")

    print("\n[ADN POLIMERASA I]")
    print("Elimina el cebador de ARN y lo sustituye por ADN.")

    print("\nCadena líder final:")
    print("5' -", abreviar(lider), "- 3'")

    # --------------------------------------------------------
    # Cadena rezagada
    # --------------------------------------------------------

    print("\n--- CADENA REZAGADA ---")
    print("Se sintetiza discontinuamente mediante fragmentos de Okazaki.")
    print("Cada fragmento se sintetiza en dirección 5' -> 3' ")
    print("La primasa coloca un cebador de ARN para iniciar cada fragmento.")

    fragmentos = []
    for i in range(0, len(adn), TAM_OKAZAKI):
        trozo_molde = adn[i:i + TAM_OKAZAKI]
        fragmento = complementaria_inversa(trozo_molde)
        fragmentos.append(fragmento)

    print("\nNúmero de fragmentos de Okazaki:", len(fragmentos))

    # Mostramos todos si la secuencia es pequeña.
    # En secuencias grandes mostramos solo los primeros.

    limite = len(fragmentos) if len(fragmentos) <= 15 else 5

    for i in range(limite):
        fragmento = fragmentos[i]
        cebador = cebador_arn(fragmento)

        print(f"\nFragmento de Okazaki {i + 1}:")
        print("  PRIMASA -> sintetiza cebador de ARN: 5'-" + cebador + "-3'")
        print("  ADN polimerasa III -> extiende el cebador: 5'-" + fragmento + "-3'")
        print("  ADN polimerasa I -> elimina el cebador y lo sustituye por ADN")

    if len(fragmentos) > limite:
        print("\n... Se omite el resto de fragmentos para facilitar la lectura.")

    print("\n[ADN LIGASA]")
    print("Une los fragmentos de Okazaki formando una cadena continua.")

    rezagada = "".join(reversed(fragmentos))

    print("\nCadena rezagada final:")
    print("5' -", abreviar(rezagada), "- 3'")

    # --------------------------------------------------------
    # Resultado
    # --------------------------------------------------------

    print("\nRESULTADO:")
    print("Se obtienen dos moléculas de ADN mediante replicación")
    print("semiconservativa: una cadena original y una nueva.")

    print("\nMolécula hija 1:")
    print("Original: 3' -", abreviar(molde), "- 5'")
    print("Nueva:    5' -", abreviar(lider), "- 3'")

    print("\nMolécula hija 2:")
    print("Original: 5' -", abreviar(adn), "- 3'")
    print("Nueva:    3' -", abreviar(rezagada[::-1]), "- 5'")

    return lider, rezagada

# ============================================================
# 5. TRANSCRIPCIÓN
# ============================================================

def transcripcion(adn):
    """Transcribe la cadena molde de ADN a ARNm."""
    titulo("PASO 2: TRANSCRIPCIÓN (ADN -> ARNm)")

    molde = complementaria(adn)

    print("\nCadena codificante:")
    print("5' -", abreviar(adn), "- 3'")

    print("\nCadena molde:")
    print("3' -", abreviar(molde), "- 5'")

    print("\n[ARN POLIMERASA]")
    print("Lee el molde 3' -> 5' y sintetiza ARNm 5' -> 3'.")

    arnm = "".join(COMPLEMENTO_ARN[base] for base in molde)

    print("\nARNm:")
    print("5' -", abreviar(arnm), "- 3'")
    print("Longitud:", len(arnm), "nucleótidos")

    return arnm

# ============================================================
# 6. TRADUCCIÓN
# ============================================================

def traduccion(arnm):
    """Traduce el ARNm desde AUG hasta un codón STOP."""
    titulo("PASO 3: TRADUCCIÓN (ARNm -> PROTEÍNA)")

    print("\nARNm:")
    print("5' -", abreviar(arnm), "- 3'")

    print("\n[RIBOSOMA]")
    print("Busca el codón de inicio AUG.")

    inicio = arnm.find(CODON_INICIO)

    if inicio == -1:
        print("\nNo existe ningún AUG. No puede comenzar la traducción.")
        return []

    print("\nAUG encontrado en la posición", inicio + 1)
    print("El ribosoma comienza la lectura de codones.\n")

    proteina = []
    resultados = []
    stop = None

    for i in range(inicio, len(arnm) - 2, 3):
        codon = arnm[i:i + 3]

        if codon in STOP_CODONES:
            stop = codon
            resultados.append((codon, "STOP"))
            break
        aminoacido = CODIGO_GENETICO.get(codon)

        if aminoacido is None:
            break

        proteina.append(aminoacido)
        resultados.append((codon, aminoacido))

    # Evitamos llenar la terminal con proteínas grandes.

    if len(resultados) <= 30:
        for numero, (codon, aa) in enumerate(resultados, 1):
            print(f"{numero:>3}. {codon} -> {aa}")
    else:
        print("La traducción contiene", len(resultados), "codones.\n")
        print("Primeros 15:")
        
        for numero, (codon, aa) in enumerate(resultados[:15], 1):
            print(f"{numero:>3}. {codon} -> {aa}")

        print("\n...")

        print("\nÚltimos 5:")
        inicio_num = len(resultados) - 4

        for numero, (codon, aa) in enumerate(
            resultados[-5:], inicio_num
        ):
            print(f"{numero:>3}. {codon} -> {aa}")

    if stop:
        print(
            f"\nCodón STOP: {stop} "
            f"({STOP_CODONES[stop]})." )

        print("Finaliza la traducción.")

    else:
        print("\nAVISO: se alcanzó el final sin encontrar un STOP")
        print("en el marco de lectura iniciado por AUG.")

    return proteina

# ============================================================
# 7. RESULTADO FINAL
# ============================================================

def resultado_final(adn, arnm, proteina):
    titulo("RESULTADO FINAL")
    print("""
FLUJO DE LA INFORMACIÓN GENÉTICA

       REPLICACIÓN
ADN ----------------> ADN

       TRANSCRIPCIÓN
ADN ----------------> ARNm

       TRADUCCIÓN
ARNm ---------------> PROTEÍNA
""")
    print("ADN:")
    print("5' -", abreviar(adn), "- 3'")

    print("\nARNm:")
    print("5' -", abreviar(arnm), "- 3'")

    print("\nPROTEÍNA:")

    if not proteina:
        print("No se ha sintetizado una proteína.")
        return

    if len(proteina) <= 25:
        print(" - ".join(proteina))

    else:
        print(
            " - ".join(proteina[:15])
            + " - ... - "
            + " - ".join(proteina[-5:]) )

    print("\nNúmero de aminoácidos:", len(proteina))

# ============================================================
# 8. PROGRAMA PRINCIPAL
# ============================================================

def main():
    titulo("SIMULADOR DEL DOGMA CENTRAL DE LA BIOLOGÍA MOLECULAR")
    print("""
El programa simula el flujo de la información genética:

ADN -> ADN -> ARNm -> PROTEÍNA
""")
    # Elegir ADN
    adn = elegir_adn()

    # Obtener la segunda cadena
    molde = complementaria(adn)

    titulo("ADN INICIAL")
    mostrar_adn(adn, molde)

    # 1. Replicación
    replicacion(adn)

    # 2. Transcripción
    arnm = transcripcion(adn)

    # 3. Traducción
    proteina = traduccion(arnm)

    # Resultado
    resultado_final(adn, arnm, proteina)

    titulo("FIN DE LA SIMULACIÓN")

# ============================================================
# EJECUCIÓN
# ============================================================

if __name__ == "__main__":
    main()