# Entropy Production Probe


# Modelo de Sondas de Büttiker para Transporte Quântico

## 📋 Índice

- [Introdução](#introdução)
- [Fundamentação Teórica](#fundamentação-teórica)
  - [Formalismo de Landauer-Büttiker](#formalismo-de-landauer-büttiker)
  - [Modelo Fenomenológico das Sondas](#modelo-fenomenológico-das-sondas)
  - [Transporte Inelástico Não-Coerente](#transporte-inelástico-não-coerente)
- [Implementação Computacional](#implementação-computacional)
  - [Sistema Modelo](#sistema-modelo)
  - [Estrutura do Código](#estrutura-do-código)
  - [Parâmetros e Configuração](#parâmetros-e-configuração)
- [Resultados e Análises](#resultados-e-análises)
- [Como Usar](#como-usar)
- [Referências](#referências)

---

## Introdução

Este código implementa o **modelo fenomenológico das sondas de Büttiker** para simular a perda de coerência de fase em transporte eletrônico quântico. O modelo é baseado no formalismo de funções de Green fora do equilíbrio (NEGF, Non-Equilibrium Green's Functions) e permite incorporar efeitos de interação elétron-fônon de forma aproximada.

### Motivação Física

Em dispositivos reais em escala nanométrica e molecular, os elétrons não se propagam de forma puramente balística. Eles interagem com as vibrações nucleares (fônons) do material, trocando energia e perdendo a coerência de fase. Essas interações:

- Causam **perda de coerência de fase** dos elétrons
- Abrem **novos canais de condutância**
- Suprimem **canais puramente elásticos**
- Dissipam energia na forma de calor


## Fundamentação Teórica

### Hamiltoniana do Sistema

O sistema completo é descrito pela Hamiltoniana:

```
H_M = H_e + H_f + H_ef
```

onde:
- **H_e**: Hamiltoniana eletrônica
- **H_f**: Hamiltoniana dos fônons (aproximação harmônica)
- **H_ef**: Hamiltoniana de interação elétron-fônon

#### Hamiltoniana Eletrônica

```
H_e = Σ_{k,σ} E_{k,σ} c†_{k,σ} c_{k,σ}
```

#### Hamiltoniana dos Fônons

```
H_f = Σ_{q,ν} ℏω_{q,ν} (b†_{q,ν} b_{q,ν} + 1/2)
```

#### Interação Elétron-Fônon

```
H_ef = Σ_{k,k',σ} Σ_{q,ν} M_{q,ν} c†_{k,σ} c_{k',σ} (b_{q,ν} + b†_{-q,ν})
```

### Formalismo de Landauer-Büttiker

O formalismo original de Landauer foi estendido por M. Büttiker para sistemas com múltiplos terminais. A corrente em um terminal p é dada por:

```
I_p = (2e²/h) Σ_q [T̄_pq V_p - T̄_qp V_q]
```

ou, no formalismo NEGF:

```
I_p = (2e/h) ∫_{-∞}^{+∞} dE Σ_{j=1}^N Tr[Γ_p G^r_M Γ†_j G^a_M] [f(E-μ_p) - f(E-μ_j)]
```

### Modelo Fenomenológico das Sondas

<img width="1678" height="970" alt="image" src="https://github.com/user-attachments/assets/8d6d8861-1c3c-4db9-b812-2479b239966b" />


**Ideia Central**: Simular o efeito do espalhamento inelástico introduzindo **sondas fictícias** (terminais) na região de espalhamento que:

1. **Acoplam-se** aos orbitais atômicos da região central
2. **Termalizam** os elétrons (perda de coerência de fase)
3. Mantêm **corrente líquida nula**: I_p = 0 para cada sonda

#### Autoenergia das Sondas

A autoenergia de interação devido às sondas é modelada como um **número imaginário puro**:

```
Σ^r_B = -iΓ/2 · S_M
```

onde:
- **Γ**: parâmetro fenomenológico de defasagem (coupling strength)
- **S_M**: matriz identidade na base da região de espalhamento

Isso produz apenas um **alargamento dos níveis** sem deslocamento de energia.

#### Função de Green Retardada

```
G^r_M = [E·S_M - H_M - Σ^r_L - Σ^r_R - Σ^r_B]^{-1}
```

### Transporte Inelástico Não-Coerente

Para satisfazer a condição I_p = 0 em cada sonda, é necessário ajustar os potenciais eletroquímicos μ_p. Isso requer resolver um sistema de **N equações não-lineares acopladas**:

```
I_p(μ_1, μ_2, ..., μ_N) = 0,  para todo p
```

#### Método de Newton para Sistemas Não-Lineares

O código implementa o método de Newton multivariável:

```
μ^(k) = μ^(k-1) - J(μ^(k-1))^{-1} · I(μ^(k-1))
```

onde **J** é a matriz Jacobiana:

```
J_ij = ∂I_i/∂μ_j
```

**Elementos da Jacobiana** (derivados analiticamente):

Para i = j:
```
J_ii = (1/2k_B T) ∫_{-∞}^{+∞} dE [T_iL + T_iR + Σ_{k≠i} T_ik] / [1 + cosh((E-μ_i)/k_B T)]
```

Para i ≠ j:
```
J_ij = -(1/2k_B T) ∫_{-∞}^{+∞} dE T_ij / [1 + cosh((E-μ_j)/k_B T)]
```

---

## Implementação Computacional

### Sistema Modelo: Cadeia Atômica 1D

![Cadeia atômica unidimensional](cap5_fig5.4.png)

O código simula uma **cadeia atômica linear** com:
- **Um nível eletrônico por sítio**
- **Acoplamento nearest-neighbor**
- **Região de espalhamento**: 1 ou mais sítios centrais
- **Eletrodos semi-infinitos**: esquerdo (L) e direito (R)

#### Vantagens deste Modelo

1. **Autoenergia analítica** para eletrodos semi-infinitos
2. **Resultados exatos** disponíveis para casos particulares
3. **Eficiência computacional**

### Autoenergia dos Eletrodos 1D

Para um eletrodo semi-infinito, a autoenergia tem forma analítica:

```python
Σ_L = (t_ML²/t_L) · exp(-ik_L)
```

onde o vetor de onda k_L é obtido da relação de dispersão:

```
E = μ_L + 2t_L cos(k_L)
⟹ k_L = arccos((E - μ_L)/(2t_L))
```

**Cuidado importante**: k_L pode ser complexo! A velocidade de grupo deve ser positiva:

```
v = ∂E/∂k_L = -2t_L sin(k_L) > 0
```

### Largura de Nível (Level Broadening)

A largura Γ associada ao acoplamento com cada eletrodo:

```python
Γ_L = -2·Im(Σ_L) = 2|sin(k_L)|
Γ_R = -2·Im(Σ_R) = 2|sin(k_R)|
```

---

## Estrutura do Código

### 1. Constantes Físicas e Parâmetros

```python
kB = 8.617333262e-5   # Constante de Boltzmann (eV/K)
T  = 10.0             # Temperatura (Kelvin)
beta = 1/(kB*T)       # Inverso da temperatura térmica
```

**Parâmetros da Tabela 5.1**:
- `eps_L, eps_M0, eps_R`: energias de sítio (eV)
- `t_L, t_R`: acoplamentos intra-eletrodo (eV)
- `t_ML, t_MR`: acoplamentos molécula-eletrodo (eV)
- `gamma_probe`: parâmetro de defasagem Γ (eV)
- `U0`: constante de Hartree (eV)

### 2. Funções Auxiliares

#### Distribuição de Fermi-Dirac

```python
def fermi(E, mu):
    """
    Função de distribuição de Fermi-Dirac
    f(E) = 1 / (1 + exp((E-μ)/(k_B T)))
    """
    x = beta*(E-mu)
    x = np.clip(x, -200, 200)  # Previne overflow
    return 1.0/(1.0 + np.exp(x))
```

#### Autoenergia de Cadeia 1D

```python
def self_energy_1D(E, eps, t, tc):
    """
    Calcula autoenergia de eletrodo semi-infinito 1D
    
    Parâmetros:
    - E: energia
    - eps: energia de sítio do eletrodo
    - t: hopping intra-eletrodo
    - tc: hopping eletrodo-molécula
    
    Retorna:
    - Σ = (tc²/t) · exp(-ik)
    """
    z = (E-eps)/(2*t)
    z = np.clip(z, -1.0, 1.0)  # Garante |z| ≤ 1
    k = np.arccos(z)
    Sigma = (tc**2/t) * np.exp(-1j*k)
    return Sigma
```

#### Largura de Nível

```python
def Gamma_from_Sigma(S):
    """
    Extrai largura Γ da autoenergia
    Γ = -2·Im(Σ)
    """
    return -2*np.imag(S)
```

### 3. Função de Green Retardada

```python
def Green(E, eps_eff):
    """
    Calcula função de Green retardada
    
    G^r = [E - ε_eff - Σ_L - Σ_R - Σ_P + iη]^{-1}
    
    onde Σ_P = -iΓ_probe/2
    """
    SL = self_energy_1D(E, eps_L, t_L, t_ML)
    SR = self_energy_1D(E, eps_R, t_R, t_MR)
    SP = -1j*gamma_probe/2  # Autoenergia da sonda
    
    Sigma = SL + SR + SP
    G = 1.0/(E - eps_eff - Sigma + 1j*eta)
    
    return G, SL, SR
```

### 4. Corrente na Sonda (Condição I_P = 0)

```python
def probe_current(muP, muL, muR, eps_eff):
    """
    Calcula corrente na sonda
    
    I_P = (2e/h) ∫ dE [T_PL(f_P-f_L) + T_PR(f_P-f_R)]
    
    onde T_ij = Γ_i·Γ_j·|G|²
    """
    Ip = 0.0
    
    for E in Egrid:
        G, SL, SR = Green(E, eps_eff)
        
        GL = Gamma_from_Sigma(SL)
        GR = Gamma_from_Sigma(SR)
        GP = gamma_probe
        
        fL = fermi(E, muL)
        fR = fermi(E, muR)
        fP = fermi(E, muP)
        
        # Coeficientes de transmissão
        TPL = GP * GL * abs(G)**2
        TPR = GP * GR * abs(G)**2
        
        # Corrente
        Ip += (TPL*(fP-fL) + TPR*(fP-fR)) * dE
    
    return Ip
```

### 5. Solver do Potencial Eletroquímico da Sonda

```python
def solve_muP(muL, muR, eps_eff):
    """
    Encontra μ_P tal que I_P(μ_P) = 0
    
    Método: Bisseção (robusto e estável)
    """
    mu_min = muR - 2.0
    mu_max = muL + 2.0
    
    for _ in range(60):  # Iterações de bisseção
        mu_mid = 0.5*(mu_min + mu_max)
        I_mid = probe_current(mu_mid, muL, muR, eps_eff)
        
        if I_mid > 0:
            mu_max = mu_mid
        else:
            mu_min = mu_mid
    
    return 0.5*(mu_min + mu_max)
```

**Observação**: Aqui foi usado bisseção ao invés do método de Newton por simplicidade, mas o método de Newton descrito na teoria também pode ser implementado para maior eficiência.

### 6. Cálculo da Ocupação (para Hartree)

```python
def occupation(muL, muR, muP, eps_eff):
    """
    Calcula ocupação eletrônica no sítio
    
    n = (1/2π) ∫ dE A(E) · [Γ_L·f_L + Γ_R·f_R + Γ_P·f_P]
    
    onde A = |G|² é a função espectral
    """
    n = 0.0
    
    for E in Egrid:
        G, SL, SR = Green(E, eps_eff)
        
        GL = Gamma_from_Sigma(SL)
        GR = Gamma_from_Sigma(SR)
        GP = gamma_probe
        
        fL = fermi(E, muL)
        fR = fermi(E, muR)
        fP = fermi(E, muP)
        
        # Função espectral ponderada
        A = abs(G)**2 * (GL*fL + GR*fR + GP*fP)
        n += A * dE/(2*pi)
    
    return n
```

### 7. Loop Principal (Sweep de Tensão)

```python
Vlist = np.linspace(0, 5.0, 60)  # Tensão de 0 a 5 V

for V in Vlist:
    # Potenciais dos eletrodos
    muL = +V/2
    muR = -V/2
    
    # 1. Resolver I_P = 0 para encontrar μ_P
    muP = solve_muP(muL, muR, eps_eff)
    
    # 2. Calcular corrente total L→R
    I = 0.0
    for E in Egrid:
        G, SL, SR = Green(E, eps_eff)
        GL = Gamma_from_Sigma(SL)
        GR = Gamma_from_Sigma(SR)
        fL = fermi(E, muL)
        fR = fermi(E, muR)
        
        TLR = GL*GR*abs(G)**2
        I += TLR*(fL-fR)*dE
    
    # 3. Atualizar energia de sítio (Hartree)
    n = occupation(muL, muR, muP, eps_eff)
    eps_eff = eps_M0 + U0*(n - n0)
```

---

## Resultados e Análises

### Efeitos do Parâmetro de Defasagem Γ

![Resultados com Γ=0.05 eV](cap5_fig5.7.png)

**Γ = 0.05 eV** (defasagem fraca):
- Alargamento moderado dos canais de condutância
- Redução da intensidade do pico de corrente
- Onset de corrente ligeiramente antecipado

**Γ = 0.11 eV** (defasagem forte):
- Alargamento significativo dos canais
- Maior supressão da corrente máxima
- Canais de condutância mais "lavados"

### Efeitos da Temperatura

![Resultados a T=300K](cap5_fig5.9.png)

**T = 10 K** vs **T = 300 K**:
- Temperatura mais alta → alargamento térmico adicional
- Distribuição de Fermi-Dirac mais suave
- Onset de corrente mais gradual
- Estruturas finas são "lavadas" pelo alargamento térmico

### Efeitos do Campo de Hartree (U₀)

**U₀ = 0.0** (sem interação elétron-elétron):
- Energia de sítio ε_M fixa
- Corrente cresce abruptamente no onset
- Comportamento mais "limpo"

**U₀ = 1.0** (com interação):
- Energia de sítio acompanha μ_L
- Ocupação fracionária gradual
- Estruturas adicionais na condutância diferencial
- Efeito espúrio do modelo (não físico)

### Efeitos do Acoplamento Molécula-Eletrodo

![Acoplamento forte](cap5_fig5.10.png)

**Acoplamento fraco** (t_ML = 0.1 eV):
- Níveis estreitos
- Onset bem definido
- Corrente surge apenas quando μ_L ≈ ε_M

**Acoplamento forte** (t_ML = 2.0 eV):
- Níveis largos
- Corrente não-nula desde V→0
- Comportamento mais "metálico"

---

## Como Usar

### Requisitos

```bash
pip install numpy matplotlib
```

### Execução Básica

```bash
python buttiker_probe_transport.py
```

### Parâmetros Ajustáveis

Modifique as variáveis no início do código:

```python
# Temperatura
T = 10.0  # ou 300.0 para temperatura ambiente

# Defasagem
gamma_probe = 0.05  # ou 0.11 para defasagem maior

# Campo de Hartree
U0 = 0.0  # ou 1.0 para incluir interação e-e

# Acoplamentos
t_ML = 0.1  # acoplamento fraco
# t_ML = 2.0  # acoplamento forte

# Grade de energia
NE = 2500  # pontos de integração
Emin, Emax = -10.0, 10.0

# Tensão
Vlist = np.linspace(0, 5.0, 60)
```

### Outputs

O código gera dois gráficos:

1. **Curva I-V**: Corrente vs. Tensão
2. **Potenciais**: μ_L, μ_R, μ_P e ε_M vs. Tensão

---

## Limitações do Modelo

⚠️ **Importante**: Este modelo fenomenológico tem limitações:

1. **Não inclui informação do espectro de fônons**
   - Γ é apenas um parâmetro global
   - Não há estrutura vibracional

2. **Não reproduz todos os efeitos experimentais**
   - Abertura de novos canais inelásticos
   - Supressão de canais elásticos
   - Estruturas finas na condutância

3. **Efeitos espúrios com U₀ ≠ 0**
   - Estruturas artificiais na condutância
   - Sonda influencia recombinação de carga

Para superar essas limitações, são necessários modelos mais sofisticados que incluam explicitamente os modos vibracionais (capítulos seguintes da tese).

---

## Referências

**Principais:**

1. **M. Büttiker**, "Four-Terminal Phase-Coherent Conductance", *Phys. Rev. Lett.* **57**, 1761 (1986)

2. **M. Büttiker**, "Symmetry of electrical conduction", *IBM J. Res. Dev.* **32**, 317 (1988)

3. **S. Datta**, *Electronic Transport in Mesoscopic Systems*, Cambridge University Press (1995)

**Implementação:**

- Capítulo 5 da tese: "Perda de Coerência de Fase"
- Seção 5.4: "Transporte Inelástico Não-Coerente"
- Seção 5.5: "Implementação"

**Método Numérico:**

- J. Stoer & R. Bulirsch, *Introduction to Numerical Analysis*, Springer (2002)

---

## Autor e Licença

Código baseado na implementação descrita no Capítulo 5 da tese sobre transporte quântico em escala molecular.

**Nota**: Este código é fornecido para fins educacionais e de pesquisa. Para uso em publicações, favor citar adequadamente.

---

## Apêndice: Caso Particular Analítico

Para **1 sítio** na região de espalhamento com **acoplamentos e energias iguais**, a transmissão L→R tem forma analítica:

```
T_LR(E) = 4·sin(k_L)·sin(k_R) / [sin(k_L) + sin(k_R)]²
```

e a ocupação da sonda (aproximação elástica):

```
f_P = [sin(k_L)·f_L + sin(k_R)·f_R] / [sin(k_L) + sin(k_R)]
```

Essas expressões foram usadas para validar o código numérico.

---

**Versão**: 1.0  
**Data**: Janeiro 2026  
**Baseado em**: Capítulo 5 - Perda de Coerência de Fase
