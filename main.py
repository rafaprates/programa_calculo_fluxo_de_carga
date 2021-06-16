##### O CÁLCULO DA DIAGONAL DA MATRIZ B ESTÁ ERRADO POR POUCO
from tabulate import tabulate
import math
import cmath
import numpy as np

# ************************************************************************* #
#                          Definição de funções                             #
# ************************************************************************* # 
def ler_arquivo(nome_arquivo):
    """
    nome_arquivo: é uma String com o endereço/nome do arquivo.
    matriz: é o arquivo lido e retornado em forma de matriz.
    """
    matriz = []
    with open(nome_arquivo, 'r') as arquivo:
        for linha in arquivo:
            linha = linha.strip('\n')
            cada_linha_do_arquivo = linha.split(' ')
            matriz.append(cada_linha_do_arquivo)
    return matriz.copy()



def extrair_coluna_da_matriz(matriz, indice_da_coluna):
    """ 
    matriz: é a matriz na qual se deseja retirar a coluna.
    indice_da_coluna: é o indice que se refere à coluna.
    """
    coluna = []
    for i in range(len(matriz)):
        elemento_da_coluna = matriz[i][indice_da_coluna]
        coluna.append(elemento_da_coluna)
    return coluna


def numerizar_vetor_string(vetor, tipo=int):
    """ 
    vetor: vetor cujos elementos são Strings
    tipo: int ou float
    vetor_numerizado: vetor mapeado para o tipo informado
    """
    vetor_numerizado = []
    vetor_numerizado = vetor.copy()
    if tipo == int:
        vetor_numerizado = map(tipo, vetor_numerizado)
    else:
        vetor_numerizado = map(tipo, vetor_numerizado)
    return list(vetor_numerizado).copy()


def contar_numero_de_elementos_unicos(vetor_a, vetor_b):
    """ 
    vetor_a: vetor de única dimensão
    vetor_b: vetor de única dimensão
    quantidade_elementos_unicos: quantidade de elementos únicos em vetor_a e vetor_b
    """
    vetor_a = vetor_a.copy()
    vetor_b = vetor_b.copy()
    quantidade_elementos_unicos = 0
    for i in range(len(vetor_a)):
        elemento = vetor_a.pop()
        if (elemento not in vetor_a) and (elemento not in vetor_b):
            quantidade_elementos_unicos += 1

    for i in range(len(vetor_b)):
        elemento = vetor_b.pop()
        if (elemento not in vetor_a) and (elemento not in vetor_b):
            quantidade_elementos_unicos += 1  

    return quantidade_elementos_unicos


def cria_matriz_vazia_mxn(M, N):
    """
    M: número de linhas
    N: número de colunas
    matriz: matriz vazia MxN
    """
    matriz = []
    for i in range(M):
        matrizAux = []
        for j in range(N):  
           matrizAux.append(0)
        matriz.append(matrizAux)
    return matriz


def vizinhos_da_barra(barra):
    vizinhos = []
    vizinhos_unicos = []
    for i in range(len(vetor_de)):
        if vetor_de[i] == barra:
            vizinhos.append(vetor_para[i])
        if vetor_para[i] == barra:
            vizinhos.append(vetor_de[i])
    for i in vizinhos:
        if i not in vizinhos_unicos:
            vizinhos_unicos.append(i)
    return vizinhos_unicos

def round_complex(x):
    return complex(round(x.real, 4),round(x.imag, 4))




# ************************************************************************* #
#                          Configurações                                    #
# ************************************************************************* # 
ARQUIVO_DELINHA = 'deLinha_livro.txt'
ARQUIVO_DEBARRA = 'deBarra_livro.txt'
MAX_ITERACAO = 1000
TOLERANCIA = 0.001
NUMERO_CASAS_DECIMAIS = 4




# ************************************************************************* #
#                          Ler Arquivos                                     #
# ************************************************************************* # 
matriz_dados_de_linha = []
matriz_dados_de_barra = []
matriz_dados_de_barra = ler_arquivo(ARQUIVO_DEBARRA)
matriz_dados_de_linha = ler_arquivo(ARQUIVO_DELINHA)

print('Dados de Linha:')
vetor_de = extrair_coluna_da_matriz(matriz_dados_de_linha, 0)
vetor_de = numerizar_vetor_string(vetor_de)
print('VETOR_DE:    ', vetor_de)

# vetor_para
vetor_para = extrair_coluna_da_matriz(matriz_dados_de_linha, 1)
vetor_para = numerizar_vetor_string(vetor_para)
print('VETOR_PARA:  ', vetor_para)

# vetor_rpu
vetor_rpu = extrair_coluna_da_matriz(matriz_dados_de_linha, 2)
vetor_rpu = numerizar_vetor_string(vetor_rpu, float)
print('VETOR_R_PU:  ', vetor_rpu)

# vetor_xpu
vetor_xpu = extrair_coluna_da_matriz(matriz_dados_de_linha, 3)
vetor_xpu = numerizar_vetor_string(vetor_xpu, float)
print('VETOR_X_PU:  ', vetor_xpu)

# vetor_bshpu
vetor_bshpu = extrair_coluna_da_matriz(matriz_dados_de_linha, 4)
vetor_bshpu = numerizar_vetor_string(vetor_bshpu, float)
print('VETOR_Bsh_PU:', vetor_bshpu)

print('\nDados de Barra:')
# vetor_P_esperado
vetor_P_esperado = extrair_coluna_da_matriz(matriz_dados_de_barra, 2)
vetor_P_esperado = numerizar_vetor_string(vetor_P_esperado, float)
print('VETOR_P_esperado:', vetor_P_esperado)

# vetor_Q_esperado
vetor_Q_esperado = extrair_coluna_da_matriz(matriz_dados_de_barra, 3)
vetor_Q_esperado = numerizar_vetor_string(vetor_Q_esperado, float)
print('VETOR_Q_esperado:', vetor_Q_esperado)

# tipo da barra
vetor_tipo_de_barra = extrair_coluna_da_matriz(matriz_dados_de_barra, 1)
print('vetor_tipo_de_barra:', vetor_tipo_de_barra)

# V da barra
vetor_V_da_barra = extrair_coluna_da_matriz(matriz_dados_de_barra, 4)
vetor_V_da_barra = numerizar_vetor_string(vetor_V_da_barra, float)
print('vetor_V_da_barra:', vetor_V_da_barra)

# Theta da barra
# V da barra
vetor_theta_da_barra = extrair_coluna_da_matriz(matriz_dados_de_barra, 5)
vetor_theta_da_barra = numerizar_vetor_string(vetor_theta_da_barra, float)
print('vetor_theta_da_barra:', vetor_theta_da_barra)
print()


barras_convergidas = []

# ************************************************************************* #
#                          Matriz Y                                         #
# ************************************************************************* # 
n_barras = contar_numero_de_elementos_unicos(vetor_de, vetor_para)
matriz_G = cria_matriz_vazia_mxn(n_barras, n_barras)
matriz_B = cria_matriz_vazia_mxn(n_barras, n_barras)


# Elementos fora da diagonal
contador_limite = range(len(vetor_rpu))
for k, m, contador in zip(vetor_de, vetor_para, contador_limite):
    r_pu = vetor_rpu[contador]
    x_pu = vetor_xpu[contador]
    
    y_pu = round_complex(-1 / complex(r_pu, x_pu))
    
    # MATRIZ G
    matriz_G[k-1][m-1] = y_pu.real
    matriz_G[m-1][k-1] = y_pu.real
    
    # MATRIZ_B
    matriz_B[k-1][m-1] = y_pu.imag
    matriz_B[m-1][k-1] = y_pu.imag


# Elementos diagonal
for lista_elementos_matriz_G, lista_elementos_matriz_B, k in zip(matriz_G, matriz_B, range(len(matriz_B))):
    # Matriz G
    matriz_G[k][k] = round(sum(lista_elementos_matriz_G) * -1, NUMERO_CASAS_DECIMAIS)
    
    # Matriz B
    matriz_B[k][k] = round(sum(lista_elementos_matriz_B) * -1, NUMERO_CASAS_DECIMAIS)
    for de, para, contador in zip(vetor_de, vetor_para, range(len(vetor_de))):
        if de-1 == k or para-1 == k:
            matriz_B[k][k] -= round(-(vetor_bshpu[contador] / 2), NUMERO_CASAS_DECIMAIS)


    

print('Matriz Y:')
print('n_barras:', n_barras)
print('\nmatriz_G:') 
print(tabulate(matriz_G))


print('\nmatriz_B:')
print(tabulate(matriz_B))

#matriz_B[0][0] = -0.9415
#matriz_B[1][1] = -0.9415
#print('\tmatriz_B:', matriz_B)
#print()

deseja_prosseguir = input("Deseja prosseguir?")
# ************************************************************************* #
#                       Calcular as potências                               #
# ************************************************************************* # 
print('Cauculo das Potências:')

# Inicializando os vetores que guardaram informações
# sobre P e Q em cada uma das barras
vetor_P_calc = []
vetor_Q_calc = []
for i in range(n_barras):
    vetor_P_calc.append(0.0)
    vetor_Q_calc.append(0.0)

# Primeiro identificar quais barra são de carga
# e quais são de referência:
vetor_barras_de_carga = []
vetor_barras_de_referencia = []
for i in range(len(vetor_tipo_de_barra)):
    if vetor_tipo_de_barra[i] == 'carga':
        vetor_barras_de_carga.append(i+1)
    else:
        vetor_barras_de_referencia.append(i+1)

print('\tvetor_barras_de_carga:',vetor_barras_de_carga)
print('\tvetor_barras_de_referencia', vetor_barras_de_referencia)

# Para as barras de carga, é necessário arbitrar um valor
# para V e outro para theta.
# Inicialmente vamos carregar esses vetores com V = 1.0 e Theta = 0.0
V_chute = 1.0
theta_chute = 0.0

# Inserir o chute inicial de V e theta no vetor
# vetor_V_da_barra e no vetor_theta_da_barra
# apenas nas posições das barras de carga
for i in range(len(vetor_tipo_de_barra)):
    if vetor_tipo_de_barra[i] == 'carga':
        vetor_V_da_barra[i] = V_chute
        vetor_theta_da_barra[i] = theta_chute

print('\nO novo vetor V e theta das barras, carregadas com o chute inicial:')
print('\tvetor_V_da_barra:', vetor_V_da_barra)
print('\tvetor_theta_da_barra', vetor_theta_da_barra)

# Agora sim eu vou calcular as potências
# --> talvez eu preciso voltar aqui
#diferenca_maior_que_tolerancia = abs(vetor_delta_P[barra-1]) <= TOLERANCIA or abs(vetor_delta_Q[barra-1]) <= TOLERANCIA
#print(diferenca_maior_que_tolerancia#)

diferenca_maior_que_tolerancia = True
#v_max = 5

v = 1
while (diferenca_maior_que_tolerancia):
    print('\n<------------------------------------------->')
    print('ITERCAO Nº:', v)
    print('<------------------------------------------->\n')
    v += 1
    for k in vetor_barras_de_carga: 
        # É necessário para loopar os vizinhos da barra e ela própria.
        # Kappa é uma lista de todas as barras vizinhas mais a própria barra
        barras_vizinhas = vizinhos_da_barra(k) 
        kappa = barras_vizinhas.copy()
        kappa.insert(0, k)

        # V_k é a tensão das barras de carga
        # V_m é a tensão das barras contidas em kappa
        # V_k = meu chute, inicialmente 1
        # theta_k = meu chute, incialmente 0
        V_k = vetor_V_da_barra[k-1]
        theta_k = vetor_theta_da_barra[k-1]

        # Antes de entrar no loop das barras vizinhas, P_calc deve ser 0!
        P_calc = 0
        Q_calc = 0
        for m in kappa:    
            V_m = vetor_V_da_barra[m-1]
            theta_m = vetor_theta_da_barra[m-1]
            G_km = matriz_G[k-1][m-1]
            B_km = matriz_B[k-1][m-1]
            theta_km = theta_k - theta_m

            P_calc += V_m * (G_km * math.cos(theta_km) + B_km * math.sin(theta_km))
            Q_calc += V_m * (G_km * math.sin(theta_km) - B_km * math.cos(theta_km))

        P_calc = V_k * P_calc
        Q_calc = V_k * Q_calc
        vetor_P_calc[k-1] = P_calc
        vetor_Q_calc[k-1] = Q_calc

    # Em posse de P_calc e Q_calc para os chutes,
    # comparar 
    vetor_delta_P = []
    vetor_delta_Q = []
    for p_esperado, p_calculado in zip(vetor_P_esperado, vetor_P_calc):
        vetor_delta_P.append(p_esperado - p_calculado)
    for q_esperado, q_calculado in zip(vetor_Q_esperado, vetor_Q_calc):
        vetor_delta_Q.append(q_esperado - q_calculado)


    print('Os resultados de Pk, Qk, delta_P e delta_Q:')
    print('\tvetor_P_calc:', vetor_P_calc)
    print('\tvetor_Q_calc:', vetor_Q_calc, '\n')
    print('\tvetor_delta_P:',vetor_delta_P)
    print('\tvetor_delta_Q:',vetor_delta_Q, '\n')

    #for delta_p, vetor_q, contador in zip(vetor_delta_P, vetor_delta_Q, range(n_barras)):
    #    vetor_barras_convergidas[contador] = True

    #for barra in vetor_barras_convergidas:
    #    if barra == False:
    #        ('\U0001F62D \U0001F62D \U0001F62D -> Não convergiu')



    erro_maximo_delta_P = abs(max(vetor_delta_P, key=abs))
    erro_maximo_delta_Q = abs(max(vetor_delta_Q, key=abs))
    
    print('erro_max_P (abs):', erro_maximo_delta_P)
    print('erro_max_Q (abs):', erro_maximo_delta_Q)


    #if 

    if erro_maximo_delta_P <= TOLERANCIA and erro_maximo_delta_Q <= TOLERANCIA:
        print('\U0001F60E \U0001F60E \U0001F60E convergiu \U0001F60E \U0001F60E \U0001F60E')
        diferenca_maior_que_tolerancia = False
    
    if v >= MAX_ITERACAO:
        diferenca_maior_que_tolerancia = False
        print('não convergiu para', v, 'iterações')

   # if ((abs(erro_maximo_delta_P) <= TOLERANCIA) and (abs(erro_maximo_delta_Q) <= TOLERANCIA)) or (MAX_ITERACAO <= v):
    #    print('\U0001F60E \U0001F60E \U0001F60E convergiu \U0001F60E \U0001F60E \U0001F60E')
    #    diferenca_maior_que_tolerancia = False
    #else:
     #   ('\U0001F62D \U0001F62D \U0001F62D -> Não convergiu')



    # ************************************************************************* #
    #                            Jacobiana                                      #
    # ************************************************************************* # 
    print('\nCalculando dV e dTheta:')
    M, N = n_barras, n_barras
    matriz_H = cria_matriz_vazia_mxn(M, N)
    matriz_N = cria_matriz_vazia_mxn(M, N)
    matriz_M = cria_matriz_vazia_mxn(M, N)
    matriz_L = cria_matriz_vazia_mxn(M, N)

    vetor_Hkk = []
    vetor_Nkk = [] 
    vetor_Mkk = []
    vetor_Lkk = []

    # Inicializando as matrizes H, N, M e L
    # em 0 para as posições das barras de referência
    for k in vetor_barras_de_referencia:
        elemento_Hkk = 0
        elemento_Nkk = 0
        elemento_Mkk = 0
        elemento_Lkk = 0

        vetor_Hkk.append(elemento_Hkk)
        vetor_Nkk.append(elemento_Nkk)
        vetor_Mkk.append(elemento_Mkk)
        vetor_Lkk.append(elemento_Lkk)

        matriz_H[k-1][k-1] = elemento_Hkk
        matriz_N[k-1][k-1] = elemento_Nkk
        matriz_M[k-1][k-1] = elemento_Mkk
        matriz_L[k-1][k-1] = elemento_Lkk

    # Inicializando as matrizes H, N, M e L
    # para as barras de carga
    for k in vetor_barras_de_carga:
        elemento_Hkk = 0
        elemento_Nkk = 0
        elemento_Mkk = 0
        elemento_Lkk = 0
        
        P_k = vetor_P_calc[k-1]
        Q_k = vetor_Q_calc[k-1]
        V_k_quadrado = vetor_V_da_barra[k-1] * vetor_V_da_barra[k-1]
        B_kk = matriz_B[k-1][k-1]
        G_kk = matriz_G[k-1][k-1]
        
        elemento_Hkk = -(B_kk * V_k_quadrado) - Q_k
        elemento_Nkk = (1/V_k) * (P_k + G_kk * V_k_quadrado)
        elemento_Mkk = -1 * (G_kk * V_k_quadrado) + P_k
        elemento_Lkk = (1/V_k) * (Q_k - B_kk * V_k_quadrado)
        
        
        vetor_Hkk.append(elemento_Hkk)
        vetor_Nkk.append(elemento_Nkk)
        vetor_Mkk.append(elemento_Mkk)
        vetor_Lkk.append(elemento_Lkk)

        matriz_H[k-1][k-1] = elemento_Hkk
        matriz_N[k-1][k-1] = elemento_Nkk
        matriz_M[k-1][k-1] = elemento_Mkk
        matriz_L[k-1][k-1] = elemento_Lkk

    print('\tvetor_H', vetor_Hkk)
    print('\tvetor_N', vetor_Nkk)
    print('\tvetor_M', vetor_Mkk)
    print('\tvetor_L', vetor_Lkk)

    matriz_coeficientes = cria_matriz_vazia_mxn(2, 2)
    matriz_termos_independentes = [0, 0]

    #for k in vetor_barras_de_referencia:
        #matriz_coeficientes[0][0] = vetor_Hkk[k-1]
        #matriz_coeficientes[0][1] = vetor_Nkk[k-1]
        #matriz_coeficientes[1][0] = vetor_Mkk[k-1]
        #matriz_coeficientes[1][1] = vetor_Lkk[k-1]

        #matriz_termos_independentes[0] = vetor_delta_P[k-1]
        #matriz_termos_independentes[1] = vetor_delta_Q[k-1]

        #incremento_theta_v = np.linalg.solve(matriz_coeficientes, matriz_termos_independentes)
        #print('inc: ', incremento_theta_v)
        #vetor_V_da_barra[k-1] += incremento_theta_v[1]
        #vetor_theta_da_barra[k-1] += incremento_theta_v[0]
    
    for k in vetor_barras_de_carga:
        print('\t\n-->Para barra ', k)
        matriz_coeficientes[0][0] = vetor_Hkk[k-1]
        matriz_coeficientes[0][1] = vetor_Nkk[k-1]
        matriz_coeficientes[1][0] = vetor_Mkk[k-1]
        matriz_coeficientes[1][1] = vetor_Lkk[k-1]

        matriz_termos_independentes[0] = vetor_delta_P[k-1]
        matriz_termos_independentes[1] = vetor_delta_Q[k-1]

        incremento_theta_v = np.linalg.solve(matriz_coeficientes, matriz_termos_independentes)
        print('inc theta: ', incremento_theta_v[0])
        vetor_V_da_barra[k-1] += incremento_theta_v[1]
        vetor_theta_da_barra[k-1] += incremento_theta_v[0]

        print('\tCoeficientes(Hkk,Nkk,Mkk,Lkk):', matriz_coeficientes)
        print('\tTermos Independentes(DeltaTheta, DeltaV):', matriz_termos_independentes)

        print('\tincremento_theta_v', incremento_theta_v)
        print('\tV das barras antes do incremento:', vetor_V_da_barra)
        print('\tV das barras apos o incremento:', vetor_V_da_barra)




print('\n\nFim <---------------------')
print('Os valores de V e Theta são:')
#print('\tvetor_V_da_barra:', vetor_V_da_barra)
#print('\tVetor_theta_da_barra:', vetor_theta_da_barra)

print('\n')
for contador, v, theta in zip(range(len(vetor_de)), vetor_V_da_barra, vetor_theta_da_barra): 
    print(f"Barra[{contador +1}] \t tensão: {round(v, NUMERO_CASAS_DECIMAIS)} \t ângulo: {theta * (180/math.pi)} graus")
