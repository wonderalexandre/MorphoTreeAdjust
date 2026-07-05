# Nota: lacunas de niveis intermediarios no Problema 1

Esta nota registra a instrumentacao adicionada ao contrato experimental do
Problema 1 para medir quando a recomposicao local nao possui um candidato
explicito para algum nivel entre `a` e `b`.

## Intervalo afetado

Para a max-tree, o evento publicado satisfaz:

```text
b > max I-(C),    a = min I-(C)
```

O intervalo de varredura local usado pela implementacao explicita e:

```text
Gamma(a,b) = {a+1, a+2, ..., b}
```

Para a min-tree dual:

```text
b < min I-(C),    a = max I-(C)
Gamma(a,b) = {a-1, a-2, ..., b}
```

Para cada `v in Gamma(a,b)`, o workspace experimental tenta observar o bucket:

```text
MergeNodes(v,-)
```

Se esse bucket esta vazio, dizemos que existe uma lacuna de candidato no
nivel `v`.

## Alvo versus intermediario

Existem dois tipos de lacuna:

```text
v = b        : falta o candidato final MergedNode(b,+)
v != b       : falta um candidato intermediario C_v^+
```

A implementacao atual ja possui uma opcao conservadora para criar o no final
`b`. Ela ainda nao cria explicitamente todos os nos intermediarios `C_v^+`.
Mesmo assim, os surveys estruturais indicam que a arvore canonica do rebuild
tem sido recuperada apos:

```text
materializacao local de C
criacao opcional do no b
varredura apenas dos buckets nao vazios
contracao canonica de arestas same-level
```

## Evidencia experimental

No marco sintetico 8x8, com validacao estrutural:

```text
eventos aceitos: 6780
interval_levels: 72806
missing_level_candidates: 11155
missing_intermediate_level_candidates: 8190
structural_mismatch: 0
```

Com a classificacao contra a arvore rebuild, apos materializar owners
parciais antes de construir o workspace diagnostico:

```text
materialization_applied: 6509
materialized_owners: 16993
materialized_nodes: 39371
workspace_candidate_levels: 61651
missing_workspace_candidate_levels: 11155
canonical_rebuild_node_levels: 63459
missing_materialized_levels: 2965
missing_intermediate_materialized_levels: 0
missing_implicit_levels: 8190
missing_intermediate_implicit_levels: 8190
NO_MATERIALIZED_INTERMEDIATE_GAP
```

No marco ICDAR 64x64 controlado:

```text
eventos aceitos: 1017
interval_levels: 113547
missing_level_candidates: 111370
missing_intermediate_level_candidates: 110359
structural_mismatch: 0
```

Com a classificacao contra a arvore rebuild, apos materializar owners
parciais antes de construir o workspace diagnostico:

```text
materialization_applied: 1017
materialized_owners: 1779
materialized_nodes: 73814
workspace_candidate_levels: 2177
missing_workspace_candidate_levels: 111370
canonical_rebuild_node_levels: 3188
missing_materialized_levels: 1011
missing_intermediate_materialized_levels: 0
missing_implicit_levels: 110359
missing_intermediate_implicit_levels: 110359
NO_MATERIALIZED_INTERMEDIATE_GAP
```

Em seguida, a busca automatica por contraexemplo foi ampliada para 100 imagens
ICDAR, dois crops 32x32 por imagem, 16 regioes conectadas por crop e
`b in {0,...,255}`:

```text
eventos aceitos: 762954
matched: 762954
missing_level_candidates: 58681270
missing_intermediate_level_candidates: 57940973
FIRST_MATERIALIZED_INTERMEDIATE_GAP: nao encontrado
NO_MATERIALIZED_INTERMEDIATE_GAP
tempo real: 444.88s
```

Portanto, lacunas intermediarias sao comuns nos eventos testados. O fato de
nao haver mismatch estrutural sugere que nem todo nivel cinza de `Gamma(a,b)`
precisa ser materializado como um no explicito da arvore final.

Mais especificamente, nesses dois marcos nenhuma lacuna intermediaria aparece
como no canonico da arvore rebuild. As lacunas que aparecem materializadas no
rebuild correspondem ao nivel alvo `b`, que a implementacao atual ja trata com
criacao opcional de no.

## Leitura teorica

O objeto `C_v^+` do paper pode ser interpretado de duas formas:

```text
objeto de prova:
  componente conexo de limiar usado para demonstrar a relacao de suporte;

no materializado:
  elemento explicito de Theta+ que aparece na arvore canonica final.
```

Os surveys favorecem a primeira leitura para muitos niveis intermediarios.
Quando dois niveis consecutivos induzem o mesmo suporte canonico, materializar
ambos criaria uma cadeia redundante que a canonicalizacao deve contrair.

A direcao teorica recomendada e formular um criterio de materializacao:

```text
criar C_v^+ somente quando o componente de limiar em v aparece como no
canonico de Theta+
```

Niveis que existem apenas no argumento de limiar, mas nao aparecem como no
canonico no rebuild, podem permanecer implicitos na prova e nao precisam virar
nos do algoritmo.

## Proxima obrigacao de prova

Precisamos provar uma equivalencia do tipo:

```text
varredura esparsa por buckets nao vazios + no b quando necessario
  -> arvore local possivelmente nao canonica, mas equivalente por reconstrucao

canonicalizacao same-level
  -> mesma arvore canonica que o rebuild de I+
```

Essa prova deve explicitar quais lacunas intermediarias sao seguras e quais
exigiriam criacao real de `C_v^+`.

O primeiro subproblema ficou mais preciso:

```text
Se v != b e MergeNodes(v,-) = vazio, entao o componente de limiar que contem C
em I+ nao aparece como no canonico de Theta+.
```

Esse enunciado ainda e uma conjectura experimental. Ate agora ele foi
observado nos marcos sintetico 8x8, ICDAR 64x64 controlado e ICDAR 100x2
region16, usando uma busca automatica por contraexemplo:

```text
print-first-materialized-intermediate-gap:
  procurar v != b tal que MergeNodes(v,-) = vazio
  e o componente de limiar em v aparece como no canonico no rebuild
```

Nenhum contraexemplo foi encontrado nesses marcos.

A formalizacao detalhada esta em
`research/docs/problem1_intermediate_gap_lemma.md`. Ela separa explicitamente:

- o componente de limiar `K_v^+(C)`;
- o objeto de prova `C_v^+`;
- o no materializado da arvore canonica;
- o corolario de completude do workspace para o caso single-`nodeCa`;
- a forma esparsa de `Hasse+` quando niveis intermediarios ficam implicitos;
- a canonicalizacao same-level como quociente canonico local.

A ponte entre essa formulacao e a mutacao concreta de
`updateTreeWithExplicitRegion` esta em
`research/docs/problem1_explicit_region_mutation_proof.md`, que agora inclui a prova do
invariante do sweep por buckets e a prova local das arestas de fronteira em
`reattachOutsideIntervalChildren`, alem do fechamento em `nodeCa` e da
absorcao de nos vazios. Ela tambem fecha a relacao entre o conjunto
operacional `S_code` e `Intervalle(C)`: o sweep visita `S_code = I_C^star`,
que contem `Intervalle(C)` e `b`; niveis extras sao apenas materializacoes de
trabalho e desaparecem por absorcao de nos vazios ou por canonicalizacao
same-level. A forma dual min-tree tambem foi explicitada nesse documento: ela
inverte as desigualdades e a ordem do sweep, mas conserva a mesma recomposicao
mutavel.
