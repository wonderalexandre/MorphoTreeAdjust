# Lema de absorcao `nodeCa -> closureRoot` no Problema 1

Este texto complementa:

```text
research/docs/problem1_general_workspace_unique_root.md
research/docs/problem1_explicit_region_mutation_proof.md
```

O objetivo e explicar por que a mutacao atual pode fechar a recomposicao em
`nodeCa`, mesmo quando a raiz teorica do workspace geral e
`closureRoot(C,b)`.

## Situacao

No diagnostico multi-anchor, definimos:

```text
closureRoot(C,b) = unico elemento de A_general(C,b).
```

Pelo teorema da raiz unica do workspace geral:

```text
closureRoot(C,b) = R_a^-(C),
```

onde `R_a^-(C)` e o representante materializado mais alto, ainda no nivel
extremo `a`, do componente conexo `K_a^-(C)` que contem `C`.

Por outro lado, o contrato operacional usado por
`updateTreeWithExplicitRegion` escolhe:

```text
max-tree: nodeCa = owner de menor altitude em O(C);
min-tree: nodeCa = owner de maior altitude em O(C).
```

Apos materializacao local, `nodeCa` e um owner vivo de algum pixel de `C` no
nivel extremo `a`. Em geral:

```text
nodeCa != closureRoot(C,b).
```

O ponto a provar e que essa diferenca nao muda a arvore canonica final quando
`nodeCa` esta no mesmo bloco same-level de `closureRoot(C,b)`.

## Lema 1: `nodeCa` pertence ao bloco same-level de `closureRoot`

Sob as hipoteses do Problema 1, apos materializacao local dos owners parciais,
`nodeCa` e descendente de `closureRoot(C,b)` e todo no no caminho ancestral de
`nodeCa` ate `closureRoot(C,b)` tem altitude `a`.

Prova. Como `nodeCa` e escolhido entre os owners de `C` com altitude extrema,
temos:

```text
alt(nodeCa) = a.
```

Tome `x in C` tal que `propernode^-(x) = nodeCa`. Como `C subset X_a(I^-)`,
o componente representado por `nodeCa` esta contido no componente conexo
`K_a^-(C)` de `X_a(I^-)` que contem `C`. Pela definicao de
`closureRoot(C,b)`, esse componente tem representante materializado superior
`R_a^-(C) = closureRoot(C,b)`.

Na arvore materializada, as divisoes internas de `K_a^-(C)` criadas para
separar proper parts internas e externas a `C` preservam a altitude `a`. Elas
podem introduzir arestas pai-filho de mesmo nivel, mas nao criam um novo
limiar entre `nodeCa` e `closureRoot(C,b)`. Logo, subindo de `nodeCa` ate o
representante superior do bloco, todos os nos encontrados continuam no nivel
`a`.

Portanto:

```text
nodeCa <=_T closureRoot(C,b)
alt(Z) = a para todo Z no caminho nodeCa -> closureRoot(C,b).
```

## Lema 2: fechar em `nodeCa` ou em `closureRoot` difere por um quociente same-level

Considere a arvore transitoria `T~+` obtida apos o sweep por buckets, antes da
canonicalizacao same-level final. Se o ultimo representante do sweep
`finalUnionNode` e fechado em `nodeCa`, entao a diferenca em relacao a fechar
diretamente em `closureRoot(C,b)` esta contida no bloco same-level:

```text
closureRoot(C,b) -> ... -> nodeCa -> finalUnionNode.
```

Todos os nos entre `closureRoot(C,b)` e `nodeCa` possuem altitude `a`.
Portanto, a contracao same-level identifica esses nos em um unico
representante de nivel `a` e reanexa `finalUnionNode` sob esse representante.

Prova. Ha dois casos operacionais.

Caso 1: `nodeCa` foi esvaziado durante a movimentacao dos pixels de `C`.
Nesse caso, `finalizeUpdateTreeAndContractRemovedNodes` substitui `nodeCa`
por `finalUnionNode` na posicao topologica de `nodeCa`, e
`absorbRemovedNodes` elimina `nodeCa` como conector vazio. Como o pai de
`nodeCa` no bloco afetado esta no mesmo caminho para `closureRoot(C,b)`, a
absorcao aproxima diretamente `finalUnionNode` da raiz same-level do bloco.
Aplicando o mesmo argumento em pos-ordem, qualquer conector vazio no caminho
ate `closureRoot(C,b)` tambem pode ser eliminado sem alterar os limiares.

Caso 2: `nodeCa` sobrevive com proper parts diretas. Nesse caso, o codigo
anexa:

```text
attachNode(nodeCa, finalUnionNode).
```

Como `nodeCa` e `closureRoot(C,b)` pertencem ao mesmo bloco de altitude `a`,
a etapa:

```text
collapseSameLevelParentChildEdges(...)
```

contrai a aresta same-level entre `nodeCa` e seu pai no caminho para
`closureRoot(C,b)`, repetindo ate que o bloco tenha um unico representante.
Quando `nodeCa` e contraido, seus proper parts e filhos, incluindo
`finalUnionNode`, passam para o representante sobrevivente do bloco. Assim, a
arvore canonicalizada e a mesma que teria sido obtida se `finalUnionNode`
tivesse sido anexado diretamente sob `closureRoot(C,b)`.

Em ambos os casos, nenhuma contracao atravessa niveis diferentes. Logo, os
componentes de limiar de `I^+` sao preservados.

## Teorema: equivalencia do fechamento operacional

Sob as hipoteses anteriores, se:

```text
nodeCa ->* closureRoot(C,b)
```

e todo o caminho esta no nivel `a`, entao fechar o sweep operacional em
`nodeCa`, seguido de absorcao de nos vazios e canonicalizacao same-level,
produz a mesma arvore canonica que fechar diretamente em
`closureRoot(C,b)`.

Esse resultado explica por que o updater atual pode continuar usando `nodeCa`
como ancora operacional local, desde que a etapa final de quociente elimine o
bloco same-level redundante.

## Evidencia experimental

O survey agora mede:

```text
nodeca_reaches_closure_root_events
nodeca_to_closure_root_same_level_path_events
closure_root_differs_but_same_level_path_events
max_nodeca_to_closure_root_path_length
```

Nos tres padroes sinteticos 6x6:

```text
accepted: 8200
matched rebuild: 8200
closure_root_differs_from_nodeca_events: 7766
nodeca_reaches_closure_root_events: 8200
nodeca_to_closure_root_same_level_path_events: 8200
closure_root_differs_but_same_level_path_events: 7766
max_nodeca_to_closure_root_path_length: 1
```

Nas duas imagens pseudoaleatorias 8x8:

```text
accepted: 6780
matched rebuild: 6780
closure_root_differs_from_nodeca_events: 6361
nodeca_reaches_closure_root_events: 6780
nodeca_to_closure_root_same_level_path_events: 6780
closure_root_differs_but_same_level_path_events: 6361
max_nodeca_to_closure_root_path_length: 1
```

Assim, nos marcos atuais, todo `nodeCa` alcanca `closureRoot(C,b)` por um
caminho same-level, e todos os casos em que `nodeCa != closureRoot(C,b)` sao
explicados por uma unica aresta same-level.

## Pendencia

Este lema resolve a diferenca entre a raiz teorica do workspace geral e a
ancora operacional escolhida pelo codigo. Ele ainda nao prova sozinho a
correcao completa para todas as regioes gerais do Problema 1. A proxima
obrigacao e mostrar que, com essa equivalencia de fechamento, os buckets
`MergeNodes_valuation(v,-)` continuam completos para todos os ramos afetados
quando `C` possui multiplos owners.
