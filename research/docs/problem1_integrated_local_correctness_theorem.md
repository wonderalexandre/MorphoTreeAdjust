# Teorema integrado de correcao local do Problema 1

Este documento integra os resultados teoricos ja registrados para o Problema
1 de atualizacao explicita por uma regiao conexa `C` e um novo nivel `b`.

Ele usa:

```text
research/docs/problem1_general_workspace_unique_root.md
research/docs/problem1_nodeca_closure_root_absorption.md
research/docs/problem1_frontier_bucket_completeness_lemma.md
research/docs/problem1_explicit_region_mutation_proof.md
```

O objetivo e enunciar uma condicao local suficiente para que a recomposicao
esparsa implementada por `updateTreeWithExplicitRegion` produza a mesma
arvore canonica que uma reconstrucao completa da imagem atualizada.

Este texto ainda nao promove o caminho experimental para producao. Ele define
o criterio teorico que deve ser verificado por asserts, diagnostics e surveys
antes dessa promocao.

## Problema

Considere uma imagem `I^-`, uma arvore de componentes canonica `T^-`, uma
regiao conexa nao vazia `C subset P` e uma imagem atualizada `I^+` tal que:

```text
I^+(x) = b, se x in C
I^+(x) = I^-(x), caso contrario.
```

Para max-tree:

```text
a = min I^-(C)
b > max I^-(C)
X_v(I) = {x in P : I(x) >= v}
Gamma(a,b) = {v : a <= v <= b}
```

Para min-tree:

```text
a = max I^-(C)
b < min I^-(C)
X_v(I) = {x in P : I(x) <= v}
Gamma(a,b) = {v : b <= v <= a}
```

O objetivo e atualizar `T^-` localmente para obter a arvore canonica de
`I^+`, sem reconstruir a arvore inteira.

## Objetos locais

Apos materializacao local dos owners cortados por `C`, definimos:

```text
O(C) = {propernode^-(x) : x in C}
```

e:

```text
W_owner^-(C,b) =
  union_{O in O(C)} ancestors_Gamma(O).
```

O workspace esparso construido pela fronteira de `C` gera os buckets:

```text
MergeNodes_valuation(v,-)
```

e denotamos:

```text
B_frontier(C,b) =
  union_v MergeNodes_valuation(v,-).
```

O workspace geral diagnosticado e:

```text
W_general^-(C,b) =
  W_owner^-(C,b) union W_sparse^-(C,b).
```

Sua raiz local unica e:

```text
closureRoot(C,b) = R_a^-(C),
```

onde `R_a^-(C)` e o representante materializado mais alto, ainda no nivel
`a`, do componente conexo `K_a^-(C)` de `X_a(I^-)` que contem `C`.

O codigo escolhe uma ancora operacional:

```text
max-tree: nodeCa = owner de menor altitude em O(C)
min-tree: nodeCa = owner de maior altitude em O(C)
```

e usa `nodeCa` para fechar a recomposicao.

## Hipoteses de liberacao

O teorema integrado vale sob as hipoteses abaixo.

```text
H1. C e nao vazio, conexo e satisfaz o Problema 1 ou seu dual min-tree.

H2. Todos os owners parciais cortados por C foram materializados localmente.

H3. A arvore materializada representa a mesma imagem I^- e so introduz
    redundancias locais removiveis por absorcao de nos vazios ou
    canonicalizacao same-level.

H4. A varredura de fronteira examina todas as adjacencias entre C e P \ C.

H5. O nivel alvo b aparece no conjunto operacional S_code, seja por bucket
    antigo, seja por criacao sintetica restrita a b.

H6. O workspace geral possui raiz unica:
    A_general(C,b) = {closureRoot(C,b)}.

H7. nodeCa alcanca closureRoot(C,b) por um caminho same-level.

H8. Nao ha gap de bucket nao explicado:
    todo X in W_owner^-(C,b) \ B_frontier(C,b) pertence ao caminho
    nodeCa -> closureRoot(C,b), ou possui Proper^-(X) != vazio e
    Proper^-(X) subset C.
```

As hipoteses H6, H7 e H8 correspondem, respectivamente, aos documentos:

```text
H6 -> research/docs/problem1_general_workspace_unique_root.md
H7 -> research/docs/problem1_nodeca_closure_root_absorption.md
H8 -> research/docs/problem1_frontier_bucket_completeness_lemma.md
```

## Teorema

Sob H1-H8, a recomposicao local que:

```text
1. materializa os owners parciais cortados por C;
2. constroi os buckets B_frontier(C,b) por varredura da fronteira;
3. cria o no alvo b quando necessario;
4. varre os buckets nao vazios em ordem dual correta;
5. move Proper^-(C) para o representante de b;
6. reanexa filhos externos e adia filhos em buckets ativos;
7. fecha a cadeia em nodeCa;
8. absorve nos vazios;
9. canonicaliza arestas same-level;
```

produz uma arvore canonica equivalente a arvore obtida por reconstrucao
completa de `I^+`.

Em termos de codigo, este e exatamente o contrato desejado para:

```text
research/include/morphoTreeAdjustResearch/DualMinMaxTreeIncrementalFilter.hpp
  updateTreeWithExplicitRegion(...)
```

quando chamado em eventos que satisfazem H1-H8.

## Prova

### 1. O workspace afetado tem uma raiz estrutural unica

Por H6, o workspace geral:

```text
W_general^-(C,b)
```

tem uma unica ancora local:

```text
closureRoot(C,b) = R_a^-(C).
```

Assim, mesmo quando `C` possui multiplos owners, a recomposicao local nao
precisa substituir uma floresta de varias raizes independentes. Todos os ramos
afetados se encontram no mesmo componente de nivel extremo `a`.

### 2. Fechar em nodeCa e equivalente a fechar em closureRoot

Por H7, `nodeCa` pertence ao mesmo bloco same-level de `closureRoot(C,b)`:

```text
nodeCa ->* closureRoot(C,b)
```

e todo no desse caminho possui altitude `a`.

Pelo lema de absorcao `nodeCa -> closureRoot`, fechar o sweep operacional em
`nodeCa`, seguido de absorcao de nos vazios e canonicalizacao same-level,
produz a mesma arvore canonica que fechar diretamente em `closureRoot(C,b)`.

Portanto, o fato de o codigo usar `nodeCa` como ancora operacional nao muda a
arvore final, desde que o caminho same-level seja realmente absorvido.

### 3. Os buckets de fronteira cobrem todos os ramos que exigem recomposicao

Por H8, se um no:

```text
X in W_owner^-(C,b)
```

nao aparece em `B_frontier(C,b)`, entao ele pertence a uma de duas classes.

Primeira classe:

```text
X esta no caminho nodeCa -> closureRoot(C,b).
```

Essa classe ja foi resolvida no passo 2: esses nos sao eliminados pelo
quociente same-level ou por absorcao de no vazio.

Segunda classe:

```text
Proper^-(X) != vazio e Proper^-(X) subset C.
```

Nesse caso, quando os pixels de `C` sao movidos para o representante de `b`,
`X` deixa de possuir proper parts diretas. Logo, `X` nao precisa de um bucket
proprio para sobreviver como componente canonico. Ele deve ser removido por
`absorbRemovedNodes` ou identificado por canonicalizacao same-level.

Consequentemente, todo no do workspace owner-interval que ainda pode carregar
proper parts externos a `C`, e que nao e apenas parte do fechamento same-level,
esta presente em algum bucket `MergeNodes_valuation(v,-)`.

### 4. O sweep esparso implementa os niveis canonicos necessarios

Pela prova de mutacao explicita, o algoritmo nao precisa visitar todos os
niveis de `Gamma(a,b)`. Ele visita:

```text
S_code = {v : MergeNodes_valuation(v,-) nao vazio} union {b}.
```

Por H5, `b in S_code`. Pela relacao entre `S_code` e `Intervalle(C)`:

```text
Intervalle(C) subset S_code.
```

Niveis intermediarios de `Gamma(a,b) \ S_code` sao implicitos: seus
componentes de limiar sao representados pelo proximo nivel materializado na
direcao dual correta. Niveis extras de `S_code \ Intervalle(C)` sao seguros,
pois seus representantes nao sobrevivem como nos canonicos novos; eles sao
absorvidos como nos vazios ou removidos por canonicalizacao same-level.

Assim, a cadeia criada pelo sweep sobre `S_code` representa a forma esparsa de
`Hasse+` para o intervalo afetado.

### 5. As arestas de fronteira sao preservadas localmente

Durante o sweep, para um no antigo `sourceNodeId` consumido pelo representante
atual `targetNodeId`, a funcao:

```text
reattachOutsideIntervalChildren(targetNodeId, sourceNodeId)
```

particiona os filhos diretos de `sourceNodeId` em:

```text
M(source) = filhos que pertencem a buckets ativos;
E(source) = filhos que nao pertencem a buckets ativos.
```

Pelo passo 3, todo filho que ainda precisa participar da recomposicao aparece
em `M(source)`. Esses filhos sao destacados e adiados para seus proprios
buckets. Os filhos em `E(source)` sao ramos externos ou ja estabilizados; ao
reanexa-los em `targetNodeId`, o codigo preserva a aresta de fronteira da
arvore atualizada.

Logo, nenhum ramo afetado e capturado prematuramente como externo, e nenhum
ramo externo e perdido pela recomposicao.

### 6. Proper, propernode e Valuation sao atualizados no lugar correto

Os representantes antigos reaproveitados preservam sua altitude. O no
sintetico de `b`, quando criado, recebe altitude `b`. As contracoes same-level
fundem apenas nos de mesma altitude. Portanto, `Valuation+` e preservada.

As operacoes:

```text
moveProperParts(...)
moveSelectedProperPartsToNode(..., C)
```

movem os proper parts dos buckets e depois movem os pixels de `C` para o
representante de `b`. Como essas operacoes atualizam `propernode` no proprio
`DynamicComponentTree`, obtemos os conjuntos `Proper+` e `propernode+`
compativeis com `I^+`.

### 7. Os conectores transitorios sao eliminados

Depois do sweep, os unicos nos que ainda podem ser redundantes sao:

```text
1. nos vazios, sem proper parts diretas;
2. nos same-level criados pela materializacao ou pelo fechamento em nodeCa.
```

Os primeiros sao removidos por:

```text
finalizeUpdateTreeAndContractRemovedNodes(...)
absorbRemovedNodes(...)
```

Os segundos sao removidos por:

```text
collapseSameLevelParentChildEdges(...)
```

Nenhuma dessas operacoes altera os componentes de limiar com pixels. Elas
apenas aplicam o quociente canonico local. Portanto, a arvore resultante e
canonica e representa `I^+`.

Combinando os passos 1-7, a arvore produzida localmente e equivalente a uma
reconstrucao completa de `I^+`.

## Criterio operacional de validacao

Antes de promover a generalizacao, o codigo deve conseguir verificar ou medir
os seguintes sinais para cada evento:

```text
accepted_by_problem1 = true
materialization_contract_valid = true
general_workspace_anchor_count = 1
nodeca_reaches_closure_root = true
nodeca_to_closure_root_same_level_path = true
owner_interval_bucket_unexplained_missing_nodes = 0
owner_interval_bucket_unexplained_gap_events = 0
target_level_present_or_created = true
```

No survey atual, os campos mais importantes sao:

```text
general_workspace_anchors
closure_root_events
nodeca_reaches_closure_root_events
nodeca_to_closure_root_same_level_path_events
owner_interval_bucket_unexplained_gap_events
owner_interval_bucket_unexplained_missing_nodes
```

O validador reutilizavel desses sinais esta em:

```text
research/include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1LocalCorrectnessValidator.hpp
```

Ele encapsula H6-H8 e o gate de nivel alvo em:

```text
validateExplicitRegionProblem1LocalCorrectnessFromAnalysis(...)
diagnoseExplicitRegionProblem1LocalCorrectness(...)
```

O updater mutavel experimental pode usar esse contrato como pre-condicao
opcional por:

```text
ExplicitRegionProblem1UpdateOptions::enforceLocalCorrectnessValidation
```

Como fronteira operacional mais conservadora, a API experimental tambem expoe
um preflight nao mutante:

```text
preflightExplicitRegionProblem1Update(...)
```

e um wrapper seguro:

```text
safeUpdateExplicitRegionProblem1(...)
```

O preflight roda contrato, materializacao diagnostica em clones, validacao
H2/H5/H6/H7/H8 e teste de capacidade do updater. O wrapper seguro so chama o
updater mutavel quando esse preflight retorna `ReadyForLocalUpdate`; caso
contrario, ele retorna o motivo estruturado e nao muta a arvore alvo.

A politica de fallback inicial tambem e experimental e opt-in:

```text
ExplicitRegionProblem1FallbackPolicy::RebuildPair
```

Quando essa politica esta habilitada, eventos que pertencem ao Problema 1
publicado, mas sao rejeitados pelo preflight local, atualizam o par max/min por
rebuild global via `RegionLevelChangeRebuildUpdater::updatePair`. Eventos
fora do contrato publicado continuam rejeitados sem mutacao.

No survey CLI, o mesmo gate e ativado por:

```text
--enforce-local-correctness
```

O resumo do survey separa as rejeicoes do gate por hipotese:

```text
local_correctness_rejected_h2
local_correctness_rejected_h5
local_correctness_rejected_h6
local_correctness_rejected_h7
local_correctness_rejected_h8
```

Esses contadores nao sao exclusivos: um mesmo evento pode violar mais de uma
hipotese.

Para localizar o primeiro caso rejeitado, use:

```text
--print-first-local-correctness-rejection
```

Um evento que viola qualquer um desses sinais deve continuar fora do contrato
de promocao. Ele deve ser tratado como contraexemplo teorico ou como falta de
instrumentacao.

Nos protocolos positivos registrados ate aqui, todos os eventos aceitos
pelo gate H1-H8 bateram tanto a reconstrucao por rebuild quanto a comparacao
estrutural canonica. O controle negativo `--no-target-node` rejeitou eventos
por H5, como esperado.

Nesse protocolo, eventos aceitos pelo Problema 1 mas rejeitados pelo preflight
local sao resolvidos por `RebuildPair` e comparados estruturalmente contra
rebuild. O stress H5 sintetico 8x8 validou 803/803 eventos aceitos, e o
stress ICDAR 10x1 crop 32x32 validou 39342/39342 eventos aceitos, ambos com
0 mismatches estruturais.

Com criacao de `b` habilitada, os primeiros crops maiores nao acionaram
fallback: ICDAR 5x1 crop 64x64 validou 10053/10053 eventos aceitos e ICDAR
2x1 crop 128x128 validou 2037/2037 eventos aceitos, ambos com 0 mismatches
estruturais. Assim, nos eventos explicitos amostrados, a questao pratica
dominante passou a ser custo do preflight, nao frequencia de fallback.

A mesma leitura apareceu nos primeiros eventos CASF reais extraidos da poda
primal.

Nesse ensaio, `C,b` sao obtidos de `rootSubtree` durante a sequencia CASF de
producao, mas o estado ainda avanca por `updateTree`. A ferramenta mede tres
camadas: um gate CASF barato, o contrato completo e o preflight diagnostico
H1-H8. O resultado atual foi:

```text
ICDAR 5x1 crop 64x64:
  755/755 eventos ready
  H6/H7/H8 = 755/755, unexplained bucket missing nodes = 0
  fast gate ~= 0.042x updateTree, 0 discordancias
  contract ~= 8.7x updateTree
  full preflight ~= 373x updateTree

ICDAR 2x1 crop 128x128:
  1282/1282 eventos ready
  H6/H7/H8 = 1282/1282, unexplained bucket missing nodes = 0
  fast gate ~= 0.043x updateTree, 0 discordancias
  contract ~= 14.9x updateTree
  full preflight ~= 1399x updateTree
```

Portanto, para esses eventos CASF, H1-H8 nao bloquearam a atualizacao local.
O bloqueio operacional e o custo da verificacao H1-H8 quando implementada
como diagnostico completo. O gate barato concordou com o preflight nesses
ensaios, e a auditoria H6-H8 nao encontrou lacunas nao explicadas. Isso ainda
nao substitui o teorema: falta provar quais hipoteses H2/H5/H6/H7/H8 sao
consequencia da estrutura CASF e quais precisam ser verificadas por asserts
baratos.

## Consequencia para implementacao

A implementacao nao precisa materializar todo `W_owner^-(C,b)`. O conjunto:

```text
B_frontier(C,b) union {nodeCa}
```

e suficiente sob H1-H8, porque os nos ausentes sao precisamente os que serao
absorvidos pelo fechamento same-level ou esvaziados quando `C` for movido para
`b`.

Isso sugere uma estrategia conservadora:

```text
1. manter updateTreeWithExplicitRegion como motor mutavel experimental;
2. usar safeUpdateExplicitRegionProblem1 como fronteira operacional;
3. ampliar surveys ate encontrar contraexemplo ou estabilizar o criterio;
4. somente entao discutir fallback e promocao do contrato geral.
```

## Pendencias abertas

O validador H2/H5/H6/H7/H8 ja existe no caminho experimental. As pendencias
agora sao de politica de uso, validacao ampliada e promocao:

```text
1. avaliar se o fallback rebuild opt-in deve permanecer como politica padrao
   experimental, ou se alguns eventos devem cair para o caminho antigo;
2. ampliar o benchmark de fallback para crops maiores e mais imagens reais;
3. revisar e consolidar
   `research/docs/problem1_casf_fast_gate_release_theorem.md` como condicao
   teorica/operacional para eventos `rootSubtree -> C,b`;
4. registrar contraexemplos minimizados quando H6, H7 ou H8 falharem;
5. decidir promocao ao CASF apenas depois da politica de fallback e dos
   custos ficarem claros.
```

Essas pendencias nao invalidam o teorema local. Elas delimitam o que ainda
precisa ser resolvido antes de transformar o backend explicito em caminho de
producao.
