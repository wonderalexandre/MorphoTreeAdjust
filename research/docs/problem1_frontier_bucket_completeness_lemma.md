# Lema de completude de fronteira dos buckets do Problema 1

Este texto transforma o diagnostico de completude dos buckets em uma
obrigacao teorica local para o Problema 1. Ele complementa:

```text
docs/problem1_general_workspace_unique_root.md
docs/problem1_nodeca_closure_root_absorption.md
```

O objetivo nao e provar que os buckets de fronteira sao iguais ao workspace
owner-interval inteiro. Essa igualdade e forte demais. O objetivo correto e
mostrar que todo no omitido dos buckets e eliminavel por uma das contracoes
locais ja previstas.

## Notacao

Trabalhamos depois da materializacao local dos owners cortados por `C`.

Para max-tree:

```text
a = min I^-(C)
b > max I^-(C)
Gamma(a,b) = {v : a <= v <= b}
```

Para min-tree:

```text
a = max I^-(C)
b < min I^-(C)
Gamma(a,b) = {v : b <= v <= a}
```

Seja:

```text
O(C) = {propernode^-(x) : x in C}
```

o conjunto dos owners diretos dos pixels de `C` na arvore materializada.
Definimos o workspace por owners como:

```text
W_owner^-(C,b) =
  union_{O in O(C)} ancestors_Gamma(O).
```

O workspace de fronteira usado pelo algoritmo e construido a partir dos
pixels `q notin C` adjacentes a algum pixel de `C`. Para cada seed
`propernode^-(q)`, a varredura sobe pelos ancestrais enquanto permanece no
lado ativo de `a`; quando o no visitado tem valuation em `Gamma(a,b)`, ele e
registrado no bucket da propria valuation:

```text
MergeNodes_valuation(v,-).
```

Denotamos por:

```text
B_frontier(C,b) =
  union_v MergeNodes_valuation(v,-)
```

o conjunto de todos os nos registrados nesses buckets.

Pelo lema de raiz unica do workspace geral, existe:

```text
closureRoot(C,b) = R_a^-(C).
```

Pelo lema de absorcao, `nodeCa` alcanca `closureRoot(C,b)` por uma cadeia
same-level nos marcos validos. Denotamos esse caminho por:

```text
P_closure(C,b) = {X : X esta no caminho nodeCa -> closureRoot(C,b)}.
```

Finalmente, `Proper^-(X)` denota os proper parts diretos de `X`, isto e, os
pixels diretamente possuidos por `X`, nao os pixels dos filhos.

## Hipoteses locais

O lema usa as seguintes hipoteses, todas alinhadas com os diagnosticos atuais:

```text
H1. C e nao vazio, conexo e aceito pelo Problema 1.

H2. Owners parciais cortados por C ja foram materializados.

H3. A arvore materializada nao possui nos vivos vazios no workspace
    owner-interval, exceto conectores same-level que serao eliminados pela
    absorcao/canonicalizacao local.

H4. A varredura de fronteira examina todas as adjacencias entre C e P \ C.

H5. nodeCa alcanca closureRoot(C,b) por um caminho same-level.
```

Sem H3, um no vivo sem proper parts diretas e fora do caminho de fechamento
seria uma terceira classe. Essa classe nao apareceu nos surveys atuais e nao
deve ser aceita silenciosamente pelo algoritmo.

## Enunciado

Lema de completude de fronteira. Sob as hipoteses do Problema 1, depois da
materializacao local dos owners cortados por `C`, seja
`X in W_owner^-(C,b)`. Se:

```text
1. X notin P_closure(C,b);
2. Proper^-(X) \ C != vazio;
```

entao:

```text
X in B_frontier(C,b).
```

Forma equivalente, mais proxima do diagnostico:

```text
Se X in W_owner^-(C,b) e X notin B_frontier(C,b),
entao X in P_closure(C,b) ou
     (Proper^-(X) != vazio e Proper^-(X) subset C).
```

Assim, um no do workspace owner-interval que esta ausente dos buckets nao e
um ramo afetado perdido. Ele e:

```text
1. um no do caminho same-level nodeCa -> closureRoot(C,b), absorvido pelo
   quociente de fechamento; ou

2. um no cujos proper parts diretos pertencem todos a C, portanto esvaziado
   quando C e movido para o representante de b.
```

## Prova

Provamos a forma direta. Seja `X in W_owner^-(C,b)` tal que:

```text
X notin P_closure(C,b)
Proper^-(X) \ C != vazio.
```

Como `X in W_owner^-(C,b)`, existe um owner `O in O(C)` tal que `X` pertence
ao caminho ancestral de `O` dentro de `Gamma(a,b)`. Logo, existe ao menos um
pixel:

```text
p in C cap support^-(X).
```

Pela segunda hipotese, existe tambem um proper part direto:

```text
r in Proper^-(X) \ C.
```

Como `support^-(X)` e conexo no limiar `Valuation^-(X)` e contem `p` e `r`,
existe um caminho em `support^-(X)` ligando `p` a `r`. Tome a primeira aresta
desse caminho que sai de `C`:

```text
c in C, q notin C, c adjacente a q.
```

O pixel `q` ainda pertence a `support^-(X)`. Fora de `C`, a imagem nao foi
modificada, entao `propernode^-(q)` e um descendente de `X` ou o proprio `X`
na arvore materializada.

A construcao dos buckets percorre exatamente esses pixels adjacentes `q`.
Como `q in support^-(X)` e `X` esta em `Gamma(a,b)`, a subida iniciada em
`propernode^-(q)` encontra `X` antes de sair do lado ativo de `a`.

Durante essa subida podem ocorrer dois tipos de nos antes de `X`:

```text
1. nos com valuation alem de b, que sao registrados como fronteira acima de b
   quando necessario;

2. nos com valuation dentro de Gamma(a,b), que sao registrados no bucket de
   sua propria valuation.
```

Quando a subida chega em `X`, temos `Valuation^-(X) in Gamma(a,b)`, pois
`X in W_owner^-(C,b)`. Portanto, a construcao chama a regra operacional
correspondente a:

```text
addMergeNode(X)
```

e registra `X` em `MergeNodes_valuation(Valuation^-(X),-)`. Assim:

```text
X in B_frontier(C,b).
```

Isso contradiz `X notin B_frontier(C,b)` e conclui a prova.

## Por que o caminho de fechamento fica fora da obrigacao

O enunciado exclui `P_closure(C,b)` porque esses nos nao precisam ser
representantes de buckets para a recomposicao ficar correta. Eles pertencem ao
bloco same-level entre a ancora operacional e a raiz teorica do workspace:

```text
nodeCa -> ... -> closureRoot(C,b).
```

Pelo lema de absorcao, fechar o sweep em `nodeCa` e depois aplicar absorcao de
nos vazios e canonicalizacao same-level produz a mesma arvore canonica que
fechar diretamente em `closureRoot(C,b)`. Portanto, a ausencia desses nos dos
buckets nao caracteriza ramo perdido; caracteriza um quociente local
redundante.

## Corolario de liberacao local

Se os dois resultados abaixo valem:

```text
1. nodeCa -> closureRoot(C,b) e um caminho same-level;
2. o lema de completude de fronteira vale para W_owner^-(C,b);
```

entao a recomposicao esparsa pode usar:

```text
B_frontier(C,b) union {nodeCa}
```

em vez de materializar todo `W_owner^-(C,b)`.

Os nos faltantes de `W_owner^-(C,b)` sao seguros porque pertencem a uma das
duas classes:

```text
closure-path:
  removidos ou identificados pelo quociente same-level;

C-only proper-part:
  esvaziados pela transferencia de C para o representante de b e removidos
  por absorbRemovedNodes.
```

Esse corolario e a forma teorica do criterio medido pelo survey:

```text
owner_interval_bucket_unexplained_gap_events = 0
owner_interval_bucket_unexplained_missing_nodes = 0
```

## Relacao com a implementacao

No codigo experimental, os objetos acima aparecem como:

```text
W_owner^-(C,b):
  collectIntervalAncestorWorkspace(...)

B_frontier(C,b):
  workspace.mergeNodesIntervalMinus

P_closure(C,b):
  isNodeOnPathToAncestor(tree, nodeCa, closureRoot, X)

Proper^-(X) subset C:
  directProperPartsAreInsideRegion(...)
```

A classificacao implementada e:

```text
X in B_frontier                       -> no coberto por bucket
X in P_closure                        -> no absorvido por fechamento
Proper^-(X) subset C                  -> no esvaziado pela mudanca de C
caso contrario                        -> gap nao explicado
```

O teste unitario atual exige que a ultima classe esteja vazia nos corpora de
survey ja usados.

## Pendencia para promocao

Este lema ainda deve ser conectado formalmente ao invariante completo do
sweep mutavel em `updateTreeWithExplicitRegion`. A proxima obrigacao e
integrar:

```text
1. raiz unica do workspace geral;
2. absorcao nodeCa -> closureRoot;
3. completude de fronteira dos buckets;
4. prova das arestas de fronteira em reattachOutsideIntervalChildren.
```

Depois dessa integracao, teremos uma condicao teorica local suficiente para
tratar regioes gerais do Problema 1 sem trocar o caminho de CASF publicado.
