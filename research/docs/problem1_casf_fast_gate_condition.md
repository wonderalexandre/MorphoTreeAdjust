# Condicao CASF para o gate barato do Problema 1

Data: 2026-06-05.

Este documento registra a primeira formalizacao operacional do gate barato
para eventos CASF reais:

```text
research/include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1CasfGate.hpp
```

O objetivo nao e substituir o teorema H1-H8. O objetivo e separar a parte que
parece ser garantida pela estrutura do evento CASF `rootSubtree -> C,b` da
parte que ainda precisa de materializacao local, criacao do nivel `b` ou
diagnostico H6-H8.

## Evento CASF

No CASF de producao, uma poda valida na arvore primal escolhe um no vivo
`rootSubtree`, diferente da raiz, e define:

```text
C = union Proper^-(X), para X no subtree(rootSubtree)
b = altitude(parent(rootSubtree))
```

Antes de podar a arvore primal, o algoritmo atualiza a arvore dual:

```text
poda na max-tree -> atualiza a min-tree
poda na min-tree -> atualiza a max-tree
```

Assim, a regiao `C` ainda esta viva na arvore primal no momento em que a
arvore dual e atualizada.

## Hipoteses CASF usadas pelo gate

O gate barato nao reconstroi workspace, nao cria clones e nao testa conexidade
por BFS. Ele assume que o chamador esta no caminho CASF e que:

```text
G1. rootSubtree e vivo, nao e a raiz e possui parent(rootSubtree).

G2. C e exatamente o suporte do subtree(rootSubtree), isto e, a uniao dos
    proper parts dos nos vivos desse subtree.

G3. C e nao vazio, sem duplicatas e conexo na adjacencia da imagem.

G4. b e a altitude de parent(rootSubtree) na arvore primal antes da poda.

G5. A arvore dual representa a mesma imagem corrente que a arvore primal.
```

Essas hipoteses sao naturais no wrapper CASF atual porque `updateTree` e
chamado antes de `pruneNode(rootSubtree)`.

## Teste barato

O gate calcula somente:

```text
minC = min I^-(C)
maxC = max I^-(C)
```

usando os owners vivos de cada pixel em `C` na arvore dual atual. Em seguida,
ele aplica a direcao do Problema 1:

```text
dual max-tree: aceitar se b > maxC
dual min-tree: aceitar se b < minC
```

No CASF, isso corresponde a:

```text
poda na min-tree:
  parent(rootSubtree) tem nivel maior que o subtree podado;
  a max-tree dual deve receber um raise local em C.

poda na max-tree:
  parent(rootSubtree) tem nivel menor que o subtree podado;
  a min-tree dual deve receber um lower local em C.
```

Esse teste e `O(|C|)` e nao aloca memoria proporcional ao dominio.

## Relacao com H1-H8

O gate cobre apenas a direcao publicada do Problema 1 sob as hipoteses CASF:

```text
H1: coberta por G1-G5 mais o teste de direcao.
```

As outras hipoteses ficam divididas em tres grupos:

```text
H2. Owners parciais cortados por C.
    Nao e garantida pelo gate. No caminho explicito atual, ela e resolvida
    pela materializacao local quando `allowLocalMaterialization=true`.

H5. Nivel alvo b presente ou criado.
    Nao e garantida pelo gate. No caminho explicito atual, ela e resolvida
    pela criacao local do no de nivel b quando
    `allowTargetNodeCreation=true`.

H6-H8. Raiz unica, fechamento same-level e completude dos buckets.
    Ainda precisam de prova especifica para eventos CASF ou de asserts
    baratos. O benchmark atual mostra 0 discordancias contra o preflight
    completo nos casos medidos, mas isso ainda e evidencia experimental.
```

Portanto, a condicao operacional candidata para promocao nao e apenas o gate:

```text
CASF fast gate
+ materializacao local habilitada
+ criacao local de b habilitada
+ prova ou asserts baratos para H6-H8
```

## Evidencia atual

Os eventos CASF reais medidos ate aqui nao expuseram falhas H6-H8. O
enunciado de liberacao CASF resume a leitura operacional dessa evidencia.

Para procurar contraexemplos H6-H8 durante a mesma sequencia CASF, use:

```text
--print-first-h6-h8-failure
```

Esse modo imprime o primeiro evento aceito pelo gate barato que falhar H6-H8,
incluindo `general_workspace_anchors`, `closure_root`, `node_ca` e lacunas de
bucket nao explicadas.

Resumo dos ensaios atuais:

```text
ICDAR 5x1 crop 64x64:
  755/755 eventos ready no preflight completo
  H6/H7/H8 = 755/755
  owner_interval_bucket_unexplained_missing_nodes = 0
  fast gate ~= 0.042x updateTree
  0 discordancias com o preflight completo

ICDAR 2x1 crop 128x128:
  1282/1282 eventos ready no preflight completo
  H6/H7/H8 = 1282/1282
  owner_interval_bucket_unexplained_missing_nodes = 0
  fast gate ~= 0.043x updateTree
  0 discordancias com o preflight completo
```

Essa evidencia sugere que o gate e barato o suficiente para hot path. O ponto
em aberto e teorico: mostrar que, para eventos CASF, H6-H8 seguem da
estrutura `rootSubtree -> C,b`, ou identificar asserts locais de custo
comparavel ao gate. A auditoria atual e favoravel: todas as lacunas de bucket
observadas foram explicadas e o workspace geral teve uma unica ancora em
todos os eventos medidos.

## Proximo passo

O proximo passo tecnico esta detalhado em:

```text
docs/problem1_casf_fast_gate_release_theorem.md
```

Esse texto organiza a prova candidata:

```text
evento CASF valido
+ fast gate
+ materializacao local
+ criacao local de b
=> H1-H8
```

Se essa prova falhar, o objetivo passa a ser encontrar os asserts baratos
minimos que separam os eventos CASF seguros dos eventos que precisam cair para
preflight completo, fallback rebuild ou `updateTree` antigo.
