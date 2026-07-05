# Teorema de liberacao CASF para o gate barato do Problema 1

Data: 2026-06-05.

Este documento formaliza a condicao que queremos usar antes de pensar em
promover o caminho explicito do Problema 1 para eventos CASF:

```text
evento CASF rootSubtree -> C,b
+ fast gate do Problema 1
+ materializacao local de C
+ criacao local do nivel b
=> hipoteses locais H1-H8
```

O objetivo e substituir o preflight H1-H8 completo no hot path por uma
condicao operacional barata, mantendo o preflight completo como diagnostico e
busca de contraexemplos.

## Referencias internas

Este texto usa os documentos:

```text
docs/problem1_casf_fast_gate_condition.md
docs/problem1_general_workspace_unique_root.md
docs/problem1_nodeca_closure_root_absorption.md
docs/problem1_frontier_bucket_completeness_lemma.md
docs/problem1_integrated_local_correctness_theorem.md
```

e o helper experimental:

```text
research/include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1CasfGate.hpp
```

## Evento CASF

Seja `Omega` o dominio da imagem. Seja `T_p^-` a arvore primal que sera
podada e `T_d^-` a arvore dual que sera atualizada. Um evento CASF valido
escolhe um no:

```text
R = rootSubtree
```

tal que `R` e vivo, nao e a raiz e possui pai `Q = parent(R)`.

Definimos:

```text
C = support^-(R)
  = union_{X in subtree(R)} Proper^-(X)

b = eta^-(Q)
```

No codigo, esse e exatamente o evento extraido antes de:

```text
updateTree(dualTree, rootSubtree)
primalTree.pruneNode(rootSubtree)
```

A ordem e importante: `C` ainda esta vivo na arvore primal no momento em que
a arvore dual e atualizada.

## Hipoteses CASF

Usaremos as seguintes hipoteses operacionais.

```text
C1. T_p^- e T_d^- representam a mesma imagem corrente I^-.

C2. R e um no vivo, diferente da raiz, com pai Q vivo.

C3. C e exatamente o suporte de R na arvore primal antes da poda.

C4. O suporte de todo no de uma arvore de componentes e conexo; logo C e
    nao vazio e conexo.

C5. O nivel de destino e b = eta^-(Q).

C6. A varredura de fronteira examina todas as adjacencias entre C e Omega \ C.

C7. O updater explicito esta autorizado a materializar owners parciais de C.

C8. O updater explicito esta autorizado a criar o no alvo de nivel b.
```

C1-C6 sao propriedades do evento CASF e do wrapper atual. C7-C8 sao
politicas do backend explicito; hoje elas correspondem a:

```text
allowLocalMaterialization = true
allowTargetNodeCreation = true
```

## Gate barato

O gate calcula, na arvore dual atual:

```text
minC = min I^-(C)
maxC = max I^-(C)
```

e aceita quando:

```text
T_d^- max-tree: b > maxC
T_d^- min-tree: b < minC
```

Esse teste e suficiente para certificar a direcao publicada do Problema 1.
No CASF, ele tambem serve como assert de sincronizacao entre a poda primal e
a atualizacao dual:

```text
poda na max-tree primal -> update da min-tree dual -> b < minC
poda na min-tree primal -> update da max-tree dual -> b > maxC
```

Se o teste falha, o evento nao deve seguir pelo caminho explicito local.

## Teorema de liberacao

Sob C1-C8, se o gate barato aceita o evento `rootSubtree -> C,b`, entao o
evento satisfaz as hipoteses H1-H8 do teorema integrado local.

Mais precisamente:

```text
H1. C e nao vazio, conexo e satisfaz o Problema 1 ou seu dual.
H2. Owners parciais cortados por C estao materializados.
H5. O nivel alvo b esta presente ou sera criado.
H6. O workspace geral possui raiz unica.
H7. nodeCa alcanca closureRoot por caminho same-level.
H8. Nao ha gap de bucket nao explicado.
```

H3 e H4 sao preservadas pelo contrato ja usado no updater: a materializacao
nao altera a imagem representada e a varredura de fronteira usa a adjacencia
completa fornecida ao CASF.

## Prova

### H1: contrato publicado

Por C2-C4, `C = support^-(R)` e nao vazio e conexo. Por C1, medir os owners
de `C` em `T_d^-` mede os niveis da imagem corrente `I^-`.

Se `T_d^-` e max-tree e o gate aceita, entao:

```text
b > max I^-(C).
```

Esse e exatamente o caso max-tree do Problema 1.

Se `T_d^-` e min-tree e o gate aceita, entao:

```text
b < min I^-(C).
```

Esse e exatamente o dual min-tree do Problema 1.

Logo H1 vale.

### H2: materializacao de owners parciais

O gate barato nao prova H2 sozinho. H2 e uma obrigacao operacional do backend
explicito.

Por C7, antes da mutacao local o updater executa a materializacao dos owners
parciais de `C`. Essa operacao separa proper parts internas e externas a `C`
sem alterar a imagem `I^-`. Depois dela, o contrato local reavaliado satisfaz:

```text
partialOwners = 0.
```

Logo H2 vale no estado materializado.

### H5: nivel alvo

O gate barato tambem nao prova H5 sozinho. H5 e resolvida por C8.

Se ja existe candidato local de valuation `b`, o sweep usa esse no. Se nao
existe, o updater explicito cria um no sintetico de nivel `b`, restrito ao
evento local. Assim:

```text
targetLevelPresentOrCreated = true.
```

Logo H5 vale.

### H6: raiz unica do workspace geral

Depois de H1 e H2, estamos exatamente nas hipoteses de:

```text
docs/problem1_general_workspace_unique_root.md
```

Como `C` e conexo, todos os owners materializados de `C` pertencem ao mesmo
componente extremo:

```text
K_a^-(C)
```

e suas cadeias ancestrais no intervalo ativo encontram o mesmo representante:

```text
closureRoot(C,b) = R_a^-(C).
```

A varredura esparsa so adiciona ramos de fronteira adjacentes a `C` que ainda
pertencem ao mesmo componente extremo. Portanto o workspace geral tem uma
unica ancora:

```text
generalWorkspaceAnchorCount = 1.
```

Logo H6 vale.

### H7: fechamento same-level

Depois de H6, `closureRoot(C,b)` e o representante materializado superior de
`K_a^-(C)` no nivel extremo `a`.

O `nodeCa` operacional e escolhido entre os owners de `C` de nivel extremo
`a`. Portanto `nodeCa` e descendente de `closureRoot(C,b)` dentro do mesmo
bloco de nivel `a`. Pela absorcao:

```text
docs/problem1_nodeca_closure_root_absorption.md
```

temos:

```text
nodeCa ->* closureRoot(C,b)
```

e todo o caminho possui valuation `a`. Fechar em `nodeCa` e depois aplicar
absorcao/canonicalizacao same-level e equivalente a fechar diretamente em
`closureRoot(C,b)`.

Logo H7 vale.

### H8: completude dos buckets de fronteira

Depois de H1-H7, aplicamos:

```text
docs/problem1_frontier_bucket_completeness_lemma.md
```

Por C6, toda adjacencia entre `C` e `Omega \ C` e examinada. Seja
`X in W_owner^-(C,b)` ausente dos buckets de fronteira. O lema mostra que
`X` pertence a uma de duas classes:

```text
1. X esta no caminho same-level nodeCa -> closureRoot(C,b);
2. Proper^-(X) != vazio e Proper^-(X) subset C.
```

A primeira classe e removida pelo quociente same-level de H7. A segunda
classe perde seus proper parts quando `C` e movido para o representante de
`b`, e e removida por absorcao de nos vazios.

Assim, nenhum no ausente dos buckets representa ramo externo perdido:

```text
owner_interval_bucket_unexplained_missing_nodes = 0.
```

Logo H8 vale.

## Condicao operacional candidata

A condicao candidata para usar o backend explicito em eventos CASF e:

```text
validPruneRoot(rootSubtree)
extract C = support(rootSubtree)
extract b = altitude(parent(rootSubtree))
checkCasfSubtreeProblem1FastGate(dualTree, C, b).readyUnderCasfAssumptions
allowLocalMaterialization = true
allowTargetNodeCreation = true
canonicalizacao same-level habilitada
absorcao de nos vazios habilitada
```

Essa condicao nao inclui o preflight H1-H8 completo. O preflight completo deve
continuar existindo como ferramenta de diagnostico, benchmark e busca de
contraexemplos.

## Asserts baratos

Os asserts baratos que podem ficar proximos ao hot path sao:

```text
A1. rootSubtree valido: vivo, nao raiz, pai vivo.
A2. C extraido do subtree e nao vazio.
A3. b igual a altitude do pai primal antes da poda.
A4. todos os pixels de C possuem owner vivo na arvore dual.
A5. gate direcional aceito: b > maxC na max-tree dual ou b < minC na min-tree dual.
A6. opcoes locais habilitadas: materializacao e criacao de b.
```

Esses asserts sao `O(|C|)` ou `O(1)` alem da propria extracao de `C`. Eles
nao substituem os asserts estruturais de depuracao; apenas delimitam o
caminho operacional candidato.

## Status

Evidencia atual:

```text
ICDAR 5x1 crop 64x64:
  755/755 eventos ready
  H6/H7/H8 = 755/755
  owner_interval_bucket_unexplained_missing_nodes = 0

ICDAR 2x1 crop 128x128:
  1282/1282 eventos ready
  H6/H7/H8 = 1282/1282
  owner_interval_bucket_unexplained_missing_nodes = 0

ICDAR 10x2 crop 64x64:
  3708/3708 eventos ready
  H6/H7/H8 = 3708/3708
  owner_interval_bucket_unexplained_missing_nodes = 0

ICDAR 5x1 crop 128x128:
  3372/3372 eventos ready
  H6/H7/H8 = 3372/3372
  owner_interval_bucket_unexplained_missing_nodes = 0
```

Com isso, a pendencia deixa de ser "medir se H6-H8 falham nos primeiros
eventos CASF" e passa a ser:

```text
1. revisar criticamente a prova acima contra o codigo real de materializacao,
   bucketizacao, absorcao e canonicalizacao;
2. ampliar a busca de contraexemplos CASF com --print-first-h6-h8-failure;
3. se nenhum contraexemplo aparecer, transformar A1-A6 em um gate operacional
   experimental antes de qualquer promocao ao CASF de producao.
```
