# Prova integrada da mutacao explicita do Problema 1

Este texto faz a ponte entre a formulacao teorica do Problema 1 e a mutacao
concreta executada por:

```text
research/include/morphoTreeAdjustResearch/DualMinMaxTreeIncrementalFilter.hpp
  updateTreeWithExplicitRegion(...)
```

O objetivo e mostrar quais operacoes do codigo implementam os objetos
`MergeNodes`, `MergedNode`, `Hasse+`, `Proper+`, `propernode+` e `Valuation+`
quando a atualizacao local recebe explicitamente uma regiao conexa `C` e um
novo nivel `b`.

Este documento nao promove o caminho experimental para o CASF de producao. Ele
organiza a prova para o contrato atual:

```text
1. C e nao vazio, conexo e ja materializado quando corta owners antigos;
2. ha um unico nodeCa;
3. o evento satisfaz o Problema 1, ou seu dual min-tree;
4. a criacao local de no e restrita ao alvo b, quando habilitada;
5. niveis intermediarios sem bucket por valuation permanecem implicitos;
6. a canonicalizacao same-level e aplicada ao final do caminho explicito.
```

## Resultados usados

Este texto usa os resultados registrados em:

```text
docs/problem1_intermediate_gap_lemma.md
docs/problem1_same_level_canonicalization_note.md
```

Em particular:

- a completude do workspace single-`nodeCa`;
- o lema de lacuna intermediaria implicita;
- a forma esparsa de `Hasse+` sobre `S_code = I_C^star`, com
  `Intervalle(C) subset S_code`;
- a canonicalizacao same-level como quociente canonico local.

## Objetos do codigo

O mapeamento principal e:

```text
properPartSetC_              -> C
altitudeCa                   -> a
b                            -> nivel alvo b
nodeCa                       -> componente antigo que ancora C no nivel a
mergeNodesByLevel_           -> buckets MergeNodes_valuation(v,-)
frontierNodesAboveB          -> fronteira acima de b
currentMergeLevel            -> nivel materializado atual do sweep
currentUnionNode             -> representante mutavel de MergedNode(currentMergeLevel,+)
previousLevelUnionNode       -> representante materializado ja processado
removedMarks_                -> nos antigos consumidos pela recomposicao
lastSameLevelEdgesCollapsed_ -> medida da canonicalizacao same-level
```

O codigo nao percorre todos os niveis de `[a,b]`. Ele percorre apenas os
niveis que aparecem em `mergeNodesByLevel_`, mais o nivel `b` quando
`allowTargetNodeCreation` cria um no sintetico. Assim, o conjunto operacional
de niveis visitados e:

```text
S_code = {v : mergeNodesByLevel_.getMergedNodes(v) nao vazio}
```

Esse conjunto nao precisa ser exatamente `Intervalle(C)`. O contrato correto e:

```text
Intervalle(C) subset S_code
```

e todo nivel extra de `S_code \ Intervalle(C)` deve ser eliminado por absorcao
de no vazio ou por canonicalizacao same-level.

## Relacao entre S_code e Intervalle(C)

O paper distingue o conjunto canonico `Intervalle(C)` de um conjunto
algoritmico maior, usualmente escrito como `I_C^star`.

Para `v in [a,b-1]`, a definicao canonica e:

```text
v in Intervalle(C)
  <=> existe X in MergeNodes_valuation(v,-)
      tal que Proper^-(X) \ C != vazio.
```

E sempre:

```text
b in Intervalle(C).
```

Ja o conjunto algoritmico visita niveis cujo bucket por valuation contem nos
com proper parts antes da remocao de `C`:

```text
v in I_C^star
  <=> existe X in MergeNodes_valuation(v,-)
      tal que Proper^-(X) != vazio,
mais o nivel b.
```

Portanto:

```text
Intervalle(C) subset I_C^star.
```

No codigo, apos a materializacao local de `C`, cada no vivo da arvore antiga
tem proper parts diretas. Como `buildMergedAndNestedCollections` insere cada
no conectado a `C` no bucket de sua propria altitude, o conjunto `S_code`
coincide com essa versao operacional de `I_C^star`, salvo pelo caso em que o
nivel `b` nao tinha candidato antigo. Quando isso acontece, o caminho coberto
pela prova exige que `createTargetLevelUnionNodeIfMissing` adicione o bucket
sintetico de `b`.

Assim, no contrato provado:

```text
S_code = I_C^star
Intervalle(C) subset S_code
b in S_code.
```

Considere agora um nivel extra:

```text
v in S_code \ Intervalle(C), v != b.
```

Pela definicao de `Intervalle(C)`, todo no `X` do bucket
`MergeNodes_valuation(v,-)` satisfaz:

```text
Proper^-(X) subset C.
```

Como o sweep processa `b` primeiro na direcao do evento, os pixels de `C` sao
movidos para o representante de `b` antes que os niveis intermediarios extras
sejam fechados. Logo, quando um nivel extra `v` e analisado, seus candidatos
ou ja estao vazios e marcados em `removedMarks_`, ou ainda podem sobreviver
temporariamente apenas como conectores estruturais sem proper part canonica
fora de `C`.

Esses niveis extras sao eliminados por uma das duas etapas ja formalizadas:

```text
1. absorcao de nos vazios:
   remove representantes sem proper parts diretas;

2. canonicalizacao same-level:
   contrai representantes redundantes que nao alteram nenhum limiar.
```

Portanto, `S_code` pode ser maior que `Intervalle(C)` sem mudar a arvore
canonica final. O papel de `S_code` e ser um conjunto de trabalho seguro; o
papel de `Intervalle(C)` e descrever os niveis que sobrevivem como nos
canonicos depois dos quocientes locais.

## Pipeline de mutacao

### 1. Entrada e ancora local

`updateTreeWithExplicitRegion` copia os pixels de entrada para
`properPartSetC_`, elimina duplicatas por `pixelsInCMarks_`, encontra o owner
de cada pixel em `dualTree->getSmallestComponent(pixelId)` e escolhe `nodeCa`
pela altitude extrema:

```text
max-tree: menor altitude entre os owners de C
min-tree: maior altitude entre os owners de C
```

Esse passo fixa `C`, `a = altitudeCa` e `nodeCa`. A prova depende da
pre-condicao de que `C` ja foi materializado quando corta owners antigos; sem
isso, `getSmallestComponent(pixelId)` pode misturar proper parts de dentro e
fora de `C` no mesmo owner.

### 2. Workspace e buckets por valuation

`buildMergedAndNestedCollections` percorre a fronteira grafica de `C`:

```text
para p em C:
  para q vizinho de p, q fora de C:
    seed = propernode^-(q)
    sobe por Hasse^- enquanto o nivel permanece no intervalo de C
```

Quando o no visitado esta no intervalo de merge, o codigo chama:

```text
mergeNodesByLevel_.addMergeNode(*dualTree, nodeSubtree)
```

Cada no entra somente no bucket da propria altitude. Portanto, esses buckets
implementam a leitura por valuation:

```text
MergeNodes_valuation(v,-)
```

Quando a subida encontra a primeira raiz de ramo acima de `b`, o codigo chama:

```text
mergeNodesByLevel_.addFrontierNodeAboveB(nodeSubtree)
```

Esse conjunto corresponde aos nos preservados que ainda precisam ser
reanexados quando o sweep atingir `b`.

### 3. Criacao opcional do alvo b

Se `allowTargetNodeCreation` esta habilitado e o bucket de `b` esta vazio,
`createTargetLevelUnionNodeIfMissing` cria:

```text
targetNode = createDetachedNode(b)
mergeNodesByLevel_.addSyntheticMergeNode(b, targetNode)
```

Esse passo implementa o caso em que `MergedNode(b,+)` precisa existir mesmo
sem candidato antigo no nivel `b`. A criacao e restrita ao alvo `b`; ela nao
cria automaticamente niveis intermediarios.

### 4. Ordem do sweep

`firstMergeLevel` compacta os buckets nao vazios. A ordem e:

```text
max-tree: niveis decrescentes
min-tree: niveis crescentes
```

No caso max-tree, isso significa processar primeiro o maior nivel materializado
e depois descer ate `a`. Se um nivel intermediario nao aparece em
`S_code`, ele e simplesmente pulado; a proxima iteracao usa diretamente o
representante materializado seguinte. Esse e o comportamento operacional da
funcao `rep(t)` definida na prova esparsa de `Hasse+`.

### 5. Fusao no bucket corrente

Em cada bucket, o primeiro no vivo nao removido vira:

```text
currentUnionNode
```

Os outros nos do mesmo bucket sao consumidos por:

```text
reattachOutsideIntervalChildren(dualTree, currentUnionNode, nodeId)
moveProperParts(currentUnionNode, nodeId)
disconnect(nodeId, release=true)
```

Interpretacao:

- `currentUnionNode` e o representante mutavel de `MergedNode(v,+)`;
- os proper parts dos nos em `MergeNodes_valuation(v,-)` sao agregados nesse
  representante;
- filhos que ainda pertencem ao intervalo afetado sao destacados para serem
  tratados por seus proprios buckets;
- filhos fora do intervalo sobem junto com o representante, preservando as
  incidencias de fronteira que dependem desse ramo.

Se todos os candidatos do bucket ja foram removidos, o codigo contrai os nos
pendentes no pai antigo e faz:

```text
currentUnionNode = previousLevelUnionNode
```

Esse e o caso operacional de nivel implicito: nao ha novo representante em
`v`, e o representante materializado anterior continua carregando o componente
de limiar.

### 6. Entrada de C no nivel b

Quando `currentMergeLevel == b`, o codigo executa:

```text
moveSelectedProperPartsToNode(dualTree, currentUnionNode, properPartSetC_)
```

Esse e o ponto em que os pixels de `C` entram no no de nivel `b`. Cada pixel e
movido de seu owner atual para `currentUnionNode` por:

```text
moveProperPart(unionNode, ownerId, pixelId)
```

Assim, `propernode+(x)` passa a apontar para o no de `b` para todo `x in C`.
Se o owner antigo fica sem proper parts, ele e marcado em `removedMarks_` e
entra em `removedNodesPendingAbsorption_`.

Ainda no nivel `b`, o codigo reanexa:

```text
frontierNodesAboveB -> currentUnionNode
```

Isso implementa a regra de `Hasse+` para ramos preservados acima de `b`:

```text
X Hasse+ MergedNode(b,+)
```

na orientacao do codigo, `currentUnionNode` e o pai do ramo preservado.

### 7. Cadeia esparsa dos novos representantes

Depois de processar o bucket corrente, se existe
`previousLevelUnionNode`, o codigo faz:

```text
attachNode(currentUnionNode, previousLevelUnionNode)
```

No caso max-tree, `previousLevelUnionNode` representa o nivel materializado
mais alto ja processado, e `currentUnionNode` representa o proximo nivel
materializado abaixo dele. Portanto, a aresta criada e a forma mutavel de:

```text
MergedNode(u,+) Hasse+ MergedNode(pred_S(u),+)
```

onde `pred_S(u)` e o nivel materializado imediatamente anterior em `S_code`.
Niveis que nao pertencem a `S_code` nao geram nos nem arestas; niveis extras
de `S_code \ Intervalle(C)` sao removidos pelos quocientes locais descritos
abaixo.

No caso min-tree, a mesma leitura vale com a ordem invertida.

### 8. Fechamento em nodeCa

`finalizeUpdateTreeAndContractRemovedNodes` resolve a ligacao entre o ultimo
representante do sweep e `nodeCa`.

Existem tres casos:

```text
1. nodeCa foi esvaziado:
   o representante final ocupa a posicao topologica de nodeCa;

2. nodeCa e a raiz:
   os ramos removidos sao normalizados e uma nova raiz compativel com a ordem
   de altitude e escolhida;

3. nodeCa sobrevive:
   o representante final e anexado como filho de nodeCa.
```

Esse passo implementa o fechamento inferior do intervalo afetado. Ele e a
parte mutavel correspondente ao caso `alpha`/ancora estavel nas proposicoes de
arestas do paper.

### 9. Absorcao de nos vazios

Depois do fechamento, `absorbRemovedNodes` elimina nos marcados que ficaram sem
proper parts. Para no nao raiz:

```text
mergedParentAndChildren(parent, removedNode)
disconnect(removedNode, release=true)
```

Para raiz removida, o codigo escolhe um novo representante de acordo com a
ordem de altitude e reanexa os demais filhos sob ele.

Essa etapa e necessaria para preservar o invariante estrutural:

```text
todo no vivo tem proper parts diretas
```

Ela nao cria novos componentes de limiar; apenas contrai estruturas que ficaram
sem proper parts apos a transferencia de `C` ou apos as fusoes dos buckets.

### 10. Canonicalizacao same-level

O caminho explicito encerra com:

```text
releaseEmptyDetachedChildlessNodes(dualTree)
collapseSameLevelParentChildEdges(dualTree)
assertAllAliveNodesHaveProperParts(dualTree)
```

`collapseSameLevelParentChildEdges` contrai qualquer aresta pai-filho com a
mesma altitude. Pelo lema de quociente canonico, essa contracao preserva todos
os conjuntos de limiar e remove apenas materializacoes redundantes.

## Invariante de sweep

Considere a lista compacta de niveis ativos construida por `firstMergeLevel`:

```text
L = (l_0, l_1, ..., l_m)
```

No caso max-tree, `L` esta em ordem decrescente. No caso min-tree, `L` esta em
ordem crescente. A prova abaixo assume que `b in L`, seja por candidato antigo,
seja pelo no sintetico criado por `createTargetLevelUnionNodeIfMissing`. Se
`b` nao esta em `L`, o evento nao esta coberto pela prova do caminho explicito,
pois `C` nao entra no representante alvo.

O sweep mantem o ponteiro:

```text
previousLevelUnionNode
```

como representante vivo do ultimo nivel materializado processado. Niveis entre
dois elementos consecutivos de `L` sao tratados como implicitos.

O invariante de inducao e:

Apos processar um nivel materializado `u`, `previousLevelUnionNode` representa
o componente recomposto associado a `u` ja com:

```text
1. todos os proper parts dos buckets anteriores a u na ordem do sweep agregados;
2. todos os filhos preservados fora do intervalo reanexados;
3. todos os ramos de niveis implicitos entre u e o proximo nivel processado
   representados pelo mesmo suporte;
4. a cadeia esparsa dos representantes ja conectada.
```

Ao entrar no proximo bucket `w = pred_S(u)`, a chamada
`attachNode(currentUnionNode, previousLevelUnionNode)` realiza a aresta:

```text
MergedNode(u,+) Hasse+ MergedNode(w,+)
```

Se nao existe bucket para um nivel `t` entre `u` e `w`, entao `t` e implicito
e nao exige mutacao propria.

### Prova do invariante

Base. Antes da primeira iteracao, nenhum bucket foi processado e
`previousLevelUnionNode = InvalidNode`. O invariante vale de forma vazia:
nenhum representante novo foi criado e nenhuma aresta da cadeia esparsa foi
introduzida.

Passo indutivo. Suponha que o invariante vale antes de processar o bucket de
nivel `u = currentMergeLevel`. O codigo zera:

```text
currentUnionNode = InvalidNode
nodesPendingRemoval_.clear()
```

e percorre todos os nos do bucket `MergeNodes_valuation(u,-)`.

#### Caso A: existe candidato vivo nao removido

O primeiro no vivo nao marcado em `removedMarks_` vira:

```text
currentUnionNode
```

Ele e destacado do pai antigo por `disconnect(..., release=false)`, porque
seu novo pai sera definido pela recomposicao local. Para cada outro no do
bucket:

```text
reattachOutsideIntervalChildren(currentUnionNode, nodeId)
moveProperParts(currentUnionNode, nodeId)
disconnect(nodeId, release=true)
```

Assim, todos os proper parts ainda vivos do bucket sao agregados em
`currentUnionNode`. Filhos que pertencem a buckets ainda nao processados sao
destacados, pois devem ser tratados no proprio nivel. Filhos fora do intervalo
afetado sobem com o representante atual.

Se `u = b`, entao:

```text
moveSelectedProperPartsToNode(currentUnionNode, C)
frontierNodesAboveB -> currentUnionNode
```

Logo, os pixels de `C` entram no representante alvo e os ramos preservados
acima de `b` passam a ter o pai prescrito pela regra esparsa de `Hasse+`.

Se ja existe `previousLevelUnionNode`, o codigo executa:

```text
attachNode(currentUnionNode, previousLevelUnionNode)
```

Pela hipotese de inducao, `previousLevelUnionNode` representa o ultimo nivel
materializado anterior na ordem do sweep. Portanto, essa chamada cria a aresta
entre representantes consecutivos de `S_code`, isto e, a forma mutavel de:

```text
MergedNode(u,+) Hasse+ MergedNode(pred_S(u),+)
```

Por fim:

```text
previousLevelUnionNode = currentUnionNode
```

Logo, o invariante vale apos o nivel `u`.

#### Caso B: nao existe candidato vivo nao removido

Nesse caso, todo no vivo encontrado no bucket ja estava marcado em
`removedMarks_`, ou entao o no ja havia sido liberado por uma iteracao
anterior. Os nos marcados sao acumulados em `nodesPendingRemoval_`.

O codigo entao executa:

```text
mergedParentAndChildren(parent(nodeId), nodeId)
disconnect(nodeId, release=true)
```

para cada no pendente ainda vivo e nao raiz. Esses nos nao possuem proper part
direta que justifique um representante novo; eles sao cascas estruturais
deixadas por movimentacoes anteriores de proper parts. Mover seus filhos para
o pai antigo preserva os componentes de limiar ja representados pela arvore
parcial, e liberar o no remove apenas uma materializacao redundante.

Depois disso o codigo apenas avanca para o proximo bucket:

```text
currentMergeLevel = nextMergeLevel()
currentUnionNode = previousLevelUnionNode
```

`previousLevelUnionNode` nao e alterado. Portanto, o nivel `u` nao cria um
novo representante e se comporta como nivel implicito ou bucket obsoleto. O
invariante permanece verdadeiro para o mesmo representante anterior.

#### Caso C: nao existe bucket para um nivel intermediario

`firstMergeLevel` e `nextMergeLevel` iteram somente sobre buckets nao vazios.
Se um nivel `t` entre dois niveis consecutivos de `L` nao aparece na lista,
entao `t notin S_code`. Para `t != b`, o lema da lacuna intermediaria
implicita mostra que `C_t^+` nao precisa ser materializado. Para `t = b`, a
pre-condicao desta prova exige que a criacao opcional de `b` tenha produzido
um bucket sintetico.

Assim, pular `t` e exatamente aplicar `rep(t)`: qualquer aresta que apontaria
para `C_t^+` aponta diretamente para o proximo representante materializado.

Conclusao. Por inducao sobre a lista compacta `L`, ao final do sweep
`previousLevelUnionNode` representa o topo da cadeia esparsa de novos
representantes, com os buckets processados agregados, as fronteiras de `b`
reanexadas e os niveis intermediarios sem representante tratados como
implicitos.

## Forma dual min-tree

O texto acima foi escrito usando a intuicao do Problema 1 max-tree. O caminho
implementado, porem, usa exatamente as mesmas operacoes mutaveis tambem para
o dual min-tree. O que muda sao apenas os predicados de ordem, a escolha de
`nodeCa` e a direcao do sweep.

No caso max-tree, o evento coberto e:

```text
a = min I^-(C)
b > max I^-(C)
X_v(I) = {x in P : I(x) >= v}
Gamma(a,b) = {a+1, ..., b}
```

No caso min-tree dual, o evento coberto e:

```text
a = max I^-(C)
b < min I^-(C)
X_v(I) = {x in P : I(x) <= v}
Gamma(a,b) = {a-1, ..., b}
```

A correspondencia operacional no codigo e:

```text
objeto                         max-tree                 min-tree
nodeCa                         owner de menor altitude   owner de maior altitude
semente valida q               altitude(q) >= a          altitude(q) <= a
subida ainda no intervalo      levelCurrent >= a         levelCurrent <= a
bucket de merge                levelCurrent <= b         levelCurrent >= b
fronteira alem de b            levelCurrent > b          levelCurrent < b
primeiro nivel do sweep        maior nivel ativo         menor nivel ativo
proximo nivel do sweep         decrementa em S_code      incrementa em S_code
condicao do laco               currentLevel > a          currentLevel < a
escolha de nova raiz           menor altitude            maior altitude
```

Assim, para min-tree, `firstMergeLevel()` inicia no menor nivel ativo. Sob a
pre-condicao desta prova, `b in S_code`; como nenhum bucket de merge e criado
abaixo de `b`, o sweep tambem processa `b` primeiro. Portanto, a entrada de
`C` em:

```text
moveSelectedProperPartsToNode(currentUnionNode, C)
```

continua acontecendo antes que os demais niveis afetados sejam fechados. Essa
e a unica propriedade de ordem usada na prova de eliminacao dos niveis extras
de `S_code \ Intervalle(C)`.

A cadeia esparsa tambem e a mesma chamada mutavel:

```text
attachNode(currentUnionNode, previousLevelUnionNode)
```

No max-tree, depois de processar `b`, o sweep visita niveis menores; logo
`currentUnionNode` tem altitude menor que `previousLevelUnionNode` e vira seu
pai. No min-tree, depois de processar `b`, o sweep visita niveis maiores; logo
`currentUnionNode` tem altitude maior que `previousLevelUnionNode` e tambem
vira seu pai, agora de acordo com a ordem da min-tree. Em ambos os casos, a
operacao cria a aresta de Hasse entre dois representantes consecutivos de
`S_code` na ordem canonica da arvore considerada.

Para os niveis implicitos, basta trocar o representante esparso:

```text
max-tree:
  rep(t) = min {u in S_code : u >= t}

min-tree:
  rep(t) = max {u in S_code : u <= t}
```

Como `b in S_code`, `rep(t)` existe no intervalo afetado em ambas as
direcoes. Se `t notin S_code`, o componente de limiar `C_t^+` e representado
pelo mesmo suporte de `C_rep(t)^+`, e nenhuma mutacao propria e necessaria no
nivel `t`.

As demais etapas sao invariantes por dualidade:

```text
1. reattachOutsideIntervalChildren
   particiona filhos por pertencerem ou nao a buckets ativos;

2. moveProperParts e moveSelectedProperPartsToNode
   atualizam Proper+ e propernode+ sem depender da ordem de altitude;

3. absorbRemovedNodes
   remove apenas nos sem proper parts diretas;

4. collapseSameLevelParentChildEdges
   contrai apenas arestas com altitudes iguais.
```

Logo, a prova integrada vale para o min-tree dual trocando `>=` por `<=`,
invertendo a ordem de `S_code` e usando a escolha dual de `nodeCa`. Nenhuma
etapa adicional de mutacao e necessaria para a forma min-tree dentro do
contrato single-`nodeCa`.

## Arestas de fronteira em reattachOutsideIntervalChildren

A funcao:

```text
reattachOutsideIntervalChildren(targetNodeId, sourceNodeId)
```

e chamada quando `sourceNodeId` sera absorvido pelo representante atual
`targetNodeId`. Seu papel e recompor apenas as arestas de fronteira locais
ligadas aos filhos diretos de `sourceNodeId`; os proper parts sao movidos por
`moveProperParts` fora dessa funcao.

Para cada filho direto `Y` de `sourceNodeId`, o codigo testa:

```text
mergeNodesByLevel_.isMergeNode(Y)
```

Esse predicado e verdadeiro exatamente quando `Y` foi colocado em algum bucket
ativo do sweep. Assim, os filhos de `sourceNodeId` sao particionados em:

```text
M(source) = {Y filho direto : Y pertence a algum bucket de merge}
E(source) = {Y filho direto : Y nao pertence a nenhum bucket de merge}
```

`M(source)` contem ramos que ainda precisam ser tratados como parte da
recomposicao local. `E(source)` contem ramos preservados em relacao ao
intervalo ativo: eles nao possuem representante proprio pendente no sweep.

### Caso sourceNodeId != targetNodeId

O codigo primeiro destaca todos os filhos marcados:

```text
if isMergeNode(Y):
    detachNode(Y)
```

Depois move todos os filhos restantes de `sourceNodeId` para `targetNodeId`:

```text
moveChildren(targetNodeId, sourceNodeId)
```

Portanto:

```text
Y in M(source)  -> Y fica temporariamente destacado;
Y in E(source)  -> Y passa a ser filho direto de targetNodeId.
```

Para `Y in E(source)`, a aresta antiga passava por `sourceNodeId`. Como
`sourceNodeId` e absorvido por `targetNodeId`, o ancestral afetado mais
proximo de `Y` na arvore atualizada passa a ser justamente `targetNodeId`.
Logo, mover `Y` para `targetNodeId` implementa a aresta de fronteira:

```text
targetNodeId Hasse+ Y
```

na orientacao da estrutura mutavel.

Para `Y in M(source)`, mover `Y` para `targetNodeId` seria prematuro. `Y`
possui um bucket proprio, logo ainda pode virar representante de
`MergedNode(w,+)` em um nivel posterior `w`. O destaque temporario impede que
um ramo interno ao intervalo seja tratado como fronteira externa. Quando o
bucket de `Y` for processado, o invariante do sweep reanexara seu
representante pela cadeia esparsa:

```text
MergedNode(w,+) Hasse+ MergedNode(pred_S(w),+).
```

Assim, a funcao nao perde a aresta: ela apenas adia a decisao ate o nivel
correto.

### Caso sourceNodeId == targetNodeId

Esse caso ocorre quando o representante escolhido no bucket ja e o proprio no
que esta sendo normalizado. O codigo nao move filhos para outro no; ele apenas
destaca os filhos marcados:

```text
if isMergeNode(Y):
    detachNode(Y)
```

Os filhos em `E(source)` permanecem sob `targetNodeId`, que ja e o
representante atual. Portanto, as arestas externas ja estao na posicao correta.
Os filhos em `M(source)` sao novamente adiados para seus buckets proprios.

### Corolario

Sob a completude do workspace single-`nodeCa`, todo filho direto que ainda
precisa participar da recomposicao aparece marcado por
`mergeNodesByLevel_.isMergeNode`. Portanto, `reattachOutsideIntervalChildren`
implementa exatamente a regra local:

```text
filho externo ao intervalo  -> reanexar ao representante atual;
filho interno ao intervalo  -> destacar e deixar para o bucket proprio.
```

Essa regra e a forma mutavel das arestas de fronteira da recomposicao esparsa
de `Hasse+`. Ela preserva os ramos estaveis e impede que ramos ainda
afetados sejam capturados pelo representante errado.

## Fechamento em nodeCa e absorcao de nos vazios

Depois do sweep, a cadeia esparsa ja foi construida, mas ainda pode haver nos
marcados por `removedMarks_` que ficaram sem proper parts. Esses nos sao
transitorios: eles nao representam uma flatzone da imagem atualizada. A etapa:

```text
finalizeUpdateTreeAndContractRemovedNodes(...)
absorbRemovedNodes(...)
```

remove esses conectores vazios e coloca o ultimo representante do sweep na
posicao topologica correta em relacao a `nodeCa`.

### Principio de quociente por no vazio

Se um no vivo `Z` esta marcado por `removedMarks_` e:

```text
Proper(Z) = vazio,
```

entao `Z` nao possui pixels cujo owner final deva ser `Z`. Seu suporte e a
uniao dos suportes dos filhos que ainda sobreviveram. Remover `Z` e reanexar
seus filhos ao representante sobrevivente mais proximo preserva a imagem
reconstruida, porque nenhum pixel tem valuation propria atribuida a `Z`.

Essa operacao e diferente da canonicalizacao same-level: aqui o motivo da
remocao nao e uma altitude repetida, mas a ausencia de proper parts diretas.
Em ambos os casos, porem, a operacao e um quociente estrutural: remove um no
transitorio sem alterar os componentes de limiar que possuem pixels.

### Caso 1: nodeCa foi removido e nao e raiz

Quando `removedMarks_.isMarked(nodeCa)` e `nodeCa` nao e raiz, o codigo faz:

```text
attachNode(parent(nodeCa), finalUnionNode)
para cada filho n de nodeCa:
    attachNode(finalUnionNode, n)
disconnect(nodeCa, release=true)
```

Como `nodeCa` foi esvaziado, ele nao tem proper parts finais. O representante
`finalUnionNode` e o topo da cadeia recomposta pelo sweep. Coloca-lo no lugar
topologico de `nodeCa` preserva o ancestral externo inferior/superior ao
intervalo afetado, dado por `parent(nodeCa)`, e mover os filhos restantes de
`nodeCa` para `finalUnionNode` preserva os ramos internos que ainda pertencem
ao componente recomposto.

Portanto, esse caso implementa a aresta de fechamento:

```text
parent(nodeCa) Hasse+ finalUnionNode
finalUnionNode Hasse+ filhos sobreviventes de nodeCa
```

na orientacao da estrutura mutavel, removendo apenas o conector vazio
`nodeCa`.

### Caso 2: nodeCa foi removido e e raiz

Quando `nodeCa` e raiz, nao ha `parent(nodeCa)` que possa receber
`finalUnionNode`. O codigo primeiro normaliza cada ramo direto por:

```text
collapseRemovedRootBranch(nodeCa, child)
```

Essa rotina promove, dentro de cada ramo, o primeiro representante que nao e
um no vazio marcado. Os demais filhos do conector vazio sao reanexados sob o
representante promovido. Assim, cada ramo exposto sob a raiz antiga passa a
comecar por um no vivo que realmente representa uma parte da arvore
atualizada.

Depois, o codigo escolhe uma nova raiz entre os filhos sobreviventes:

```text
max-tree: menor altitude
min-tree: maior altitude
```

e reanexa os demais filhos sob essa raiz. Essa escolha preserva a ordem global
da arvore de componentes: a raiz deve ser o componente total no nivel extremo.
Por fim, `nodeCa` e liberado. Como `nodeCa` estava vazio, a operacao remove
apenas a raiz transitoria e preserva os suportes dos ramos sobreviventes.

### Caso 3: nodeCa sobrevive

Se `nodeCa` nao foi marcado e `finalUnionNode != nodeCa`, o codigo executa:

```text
attachNode(nodeCa, finalUnionNode)
```

Nesse caso, `nodeCa` ainda possui proper parts e continua sendo o componente
estavel que ancora a base do intervalo afetado. O representante final do sweep
deve aparecer como filho de `nodeCa`, exatamente como o primeiro nivel
materializado do intervalo a partir de `a` na forma esparsa de `Hasse+`.

### Absorcao pos-ordem dos nos vazios

Depois do fechamento, `absorbRemovedNodes` percorre os nos marcados em
pos-ordem. A pos-ordem e importante porque um no vazio pode conter filhos que
tambem estao marcados; os filhos devem ser normalizados antes do pai.

Para no nao raiz, a operacao e:

```text
mergedParentAndChildren(parent, removedNode)
disconnect(removedNode, release=true)
```

Como `removedNode` nao tem proper parts, mover filhos e proper parts para o
pai nao muda owner final de nenhum pixel de `removedNode`: nao ha pixel direto
em `removedNode`. A operacao apenas elimina um nivel artificial da cadeia.
Se o pai tambem fica sem proper parts, ele e marcado para absorcao posterior,
propagando o quociente ate encontrar um representante vivo.

Para raiz removida, `absorbRemovedRootNode` escolhe entre os filhos
sobreviventes o representante compativel com a ordem extrema de altitude e
reanexa os demais filhos sob ele. Isso e o mesmo argumento do caso em que
`nodeCa` era raiz: uma raiz vazia nao representa flatzone; a nova raiz deve ser
o componente total sobrevivente.

### Corolario

Sob o invariante de sweep e a prova de reanexacao de fronteira,
`finalizeUpdateTreeAndContractRemovedNodes` seguido de `absorbRemovedNodes`
preserva a imagem reconstruida e os componentes de limiar ja obtidos pelo
sweep. A etapa remove apenas conectores vazios e reposiciona
`finalUnionNode` no ponto de fechamento do intervalo afetado.

Assim, antes da canonicalizacao same-level, a arvore resultante e
equivalente por reconstrucao a arvore esparsa de `Hasse+`, podendo ainda
conter somente redundancias de mesmo nivel que serao removidas por
`collapseSameLevelParentChildEdges`.

## Invariantes de Proper, propernode e Valuation

O codigo preserva `Valuation+` porque:

```text
1. no reaproveitado no bucket v ja tem altitude v;
2. no sintetico criado para b recebe altitude b;
3. contracao same-level so funde nos de mesma altitude.
```

O codigo constroi `Proper+` em duas etapas:

```text
1. moveProperParts(currentUnionNode, nodeId)
   agrega os proper parts dos nos antigos do bucket;

2. moveSelectedProperPartsToNode(..., C)
   move todos os pixels de C para o representante de b.
```

Assim, para `v != b`, qualquer pixel de `C` que tenha sido carregado
temporariamente por um representante inferior e retirado quando o nivel `b` e
processado. Para `v = b`, `C` passa a fazer parte dos proper parts do
representante alvo.

O codigo constroi `propernode+` porque as operacoes `moveProperPart` e
`moveProperParts` atualizam o owner dos pixels movidos na propria
`DynamicComponentTree`.

## Estado da prova e obrigacoes abertas

Este documento fecha o invariante do sweep por buckets, a prova local de
`reattachOutsideIntervalChildren` e o quociente por nos vazios no fechamento
de `nodeCa`, no contrato single-`nodeCa`, incluindo o caso em que
`currentUnionNode == InvalidNode`.

O uso desses argumentos dentro do contrato geral do Problema 1, com multiplos
owners de `C`, esta organizado em:

```text
docs/problem1_integrated_local_correctness_theorem.md
```

A relacao entre `S_code` e `Intervalle(C)` tambem fica fechada sob o contrato
atual: o codigo varre `S_code = I_C^star`, nao necessariamente
`Intervalle(C)`, e temos:

```text
Intervalle(C) subset S_code
b in S_code
```

Os niveis extras de `S_code \ Intervalle(C)` nao representam novos nos
canonicos. Eles sao absorvidos como nos vazios ou removidos pela
canonicalizacao same-level.

Com a secao dual min-tree acima, a prova integrada esta fechada para o escopo
single-`nodeCa` atualmente implementado.

Permanece fora do escopo atual: casos de multiplos `nodeCa` ou decomposicao
multi-anchor. Eles pertencem ao plano de regiao geral iniciado em
`docs/problem1_general_region_multi_anchor.md`.

## Leitura operacional

Com a teoria atual, o caminho correto nao e criar todos os `C_v^+` de
`Gamma(a,b)`. O caminho correto e:

```text
1. materializar C quando necessario;
2. construir buckets por valuation a partir da fronteira;
3. criar b se o alvo nao tem candidato;
4. varrer apenas buckets materializados;
5. reanexar fronteiras e cadeia esparsa;
6. absorver nos vazios;
7. canonicalizar same-level.
```

Se o oracle encontrar um caso em que essa sequencia reconstroi a imagem mas
nao bate a arvore canonica apos same-level, o contraexemplo deve ser
classificado no estudo multi-anchor ou em uma nova obrigacao do contrato
single-`nodeCa`, conforme sua estrutura.
