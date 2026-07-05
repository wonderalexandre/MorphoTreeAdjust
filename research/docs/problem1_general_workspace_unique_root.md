# Teorema da raiz unica do workspace geral no Problema 1

Este texto formaliza a observacao experimental registrada em
`docs/problem1_general_region_multi_anchor.md`: mesmo quando a regiao conexa
`C` possui multiplos owners depois da materializacao local, o workspace geral
medido pelo diagnostico possui uma unica ancora local.

O objetivo nao e provar ainda a correcao completa da mutacao local para
regioes gerais. O objetivo e fechar uma questao anterior: decidir se o caso
geral exige uma recomposicao por floresta de varias ancoras ou se ainda existe
uma raiz estrutural unica do workspace afetado.

## Notacao

Para uma max-tree, usamos:

```text
X_v(I) = {x in P : I(x) >= v}.
```

Para uma min-tree, usamos:

```text
X_v(I) = {x in P : I(x) <= v}.
```

O evento do Problema 1 e:

```text
max-tree:
  a = min I^-(C)
  b > max I^-(C)

min-tree:
  a = max I^-(C)
  b < min I^-(C)
```

Defina o intervalo ativo:

```text
Gamma(a,b) = {v : a <= v <= b}, para max-tree;
Gamma(a,b) = {v : b <= v <= a}, para min-tree.
```

Depois da materializacao local dos owners cortados por `C`, seja:

```text
O(C) = {propernode^-(x) : x in C}
```

o conjunto dos owners diretos de `C` na arvore materializada. A materializacao
nao altera a imagem `I^-`; ela apenas torna partes internas/externas de owners
parciais explicitamente representadas como nos locais, possivelmente com
arestas pai-filho de mesmo nivel que serao removidas pela canonicalizacao.

Seja `K_a^-(C)` o componente conexo de `X_a(I^-)` que contem `C`. Como `C` e
conexo e todo pixel de `C` pertence a `X_a(I^-)`, esse componente existe e e
unico.

Na arvore materializada, denotamos por `R_a^-(C)` o no vivo mais alto, ainda
no nivel `a`, que representa esse componente no bloco local de
`K_a^-(C)`. Na arvore canonica original, `R_a^-(C)` coincide com o no
canonico de `K_a^-(C)`. Depois da materializacao, podem existir filhos de
mesmo nivel abaixo desse representante, mas eles nao substituem a raiz local
do bloco.

No codigo, o workspace por owners e:

```text
W_owner^-(C,b) =
  union_{O in O(C)} ancestors_Gamma(O),
```

onde `ancestors_Gamma(O)` sobe pela relacao pai-filho da arvore materializada
enquanto a altitude do no pertence a `Gamma(a,b)`.

O workspace esparso usado pelo updater atual e denotado por:

```text
W_sparse^-(C,b).
```

O diagnostico mais geral usa:

```text
W_general^-(C,b) =
  W_owner^-(C,b) union W_sparse^-(C,b).
```

As ancoras locais sao as raizes induzidas do workspace:

```text
A_general(C,b) =
  {X in W_general^-(C,b) : parent^-(X) notin W_general^-(C,b)}.
```

## Lema 1: todos os owners sobem ate `R_a^-(C)` dentro do intervalo ativo

Para todo `O in O(C)`, o caminho ancestral de `O` contem `R_a^-(C)` antes de
sair de `Gamma(a,b)`.

Prova para max-tree. Tome `x in C` tal que `propernode^-(x) = O`. Pela
definicao de `a`, temos `I^-(x) >= a`, entao `x in X_a(I^-)`. Como `C` e
conexo e `C subset X_a(I^-)`, todos os pixels de `C` pertencem ao mesmo
componente conexo `K_a^-(C)` de `X_a(I^-)`.

O no `O` representa um componente de algum nivel `h >= a` que contem `x`.
Pela propriedade de inclusao da max-tree, esse componente esta contido no
componente de nivel `a` que contem `x`, isto e, em `K_a^-(C)`. Portanto, no
caminho de `O` para a raiz da arvore, encontramos o representante vivo
`R_a^-(C)`. As altitudes ao longo desse trecho pertencem ao intervalo
`[a,h]`, e `h <= max I^-(C) < b`. Logo, todo esse trecho pertence a
`Gamma(a,b)`.

O argumento para min-tree e dual. Nesse caso `a = max I^-(C)`, todo pixel de
`C` pertence a `X_a(I^-) = {x : I^-(x) <= a}`, e cada owner de nivel
`h <= a` esta contido no componente de nivel `a` que contem `C`. Como
`b < min I^-(C) <= h`, o trecho de `O` ate `R_a^-(C)` pertence a
`Gamma(a,b)`.

Portanto:

```text
R_a^-(C) in W_owner^-(C,b)
```

e todo owner de `C` esta conectado a esse mesmo no dentro de
`W_owner^-(C,b)`.

## Corolario 1: `W_owner^-(C,b)` possui uma unica ancora

Pelo Lema 1, cada cadeia ancestral iniciada em um owner de `C` encontra o
mesmo no `R_a^-(C)` antes de sair do intervalo ativo. Como o workspace
`W_owner^-(C,b)` e a uniao dessas cadeias, ele e conectado como subgrafo
induzido da arvore materializada.

Nenhum ancestral estrito de `R_a^-(C)` pertence a `W_owner^-(C,b)`. Pela
definicao de `R_a^-(C)`, ele e o representante mais alto do bloco local de
nivel `a`. Para a max-tree, qualquer ancestral acima desse representante tem
altitude menor que `a`; para a min-tree, qualquer ancestral acima desse
representante tem altitude maior que `a`. Em ambos os casos, esse ancestral
esta fora de `Gamma(a,b)`.

Assim:

```text
|A_owner(C,b)| = 1.
```

## Lema 2: o workspace esparso nao introduz uma nova ancora

Todo no de `W_sparse^-(C,b)` tambem esta abaixo de `R_a^-(C)` ou coincide com
ele, quando restrito ao intervalo ativo.

No codigo, `W_sparse^-(C,b)` e construido a partir de sementes fora de `C`
que sao adjacentes a algum pixel de `C`. Considere uma semente `q notin C`
adjacente a `p in C`.

Para max-tree, a semente so e aceita pelo workspace se seu owner tem altitude
`h >= a`; logo `q in X_a(I^-)`. Como `p in C subset X_a(I^-)` e `p ~ q`, os
pixels `p` e `q` pertencem ao mesmo componente conexo de `X_a(I^-)`, isto e,
ao componente `K_a^-(C)`. Portanto, o owner de `q` e todos os seus ancestrais
coletados dentro de `Gamma(a,b)` estao no ramo que encontra `R_a^-(C)`.

Para min-tree, a semente so e aceita quando seu owner tem altitude `h <= a`;
logo `q in X_a(I^-)`. O mesmo argumento de adjacencia mostra que `q` pertence
ao componente `K_a^-(C)`, e a cadeia coletada encontra `R_a^-(C)`
antes de sair de `Gamma(a,b)`.

Consequentemente, ao adicionar `W_sparse^-(C,b)` a `W_owner^-(C,b)`, nenhum
ramo novo fica desconectado da raiz unica ja presente em `W_owner^-(C,b)`.

## Teorema: raiz unica do workspace geral

Sob as hipoteses do Problema 1, com `C` conexo e apos materializacao local dos
owners parciais de `C`, vale:

```text
|A_general(C,b)| = 1.
```

Mais precisamente, a unica ancora de `W_general^-(C,b)` e `R_a^-(C)`, o
representante materializado do componente `K_a^-(C)` no nivel extremo `a`.

Prova. Pelo Corolario 1, `W_owner^-(C,b)` possui uma unica ancora e contem o
representante `R_a^-(C)`. Pelo Lema 2, todo no adicionado por
`W_sparse^-(C,b)` esta em uma cadeia que encontra esse mesmo representante
dentro de `Gamma(a,b)`. Portanto a uniao
`W_general^-(C,b) = W_owner^-(C,b) union W_sparse^-(C,b)` e conectada e tem
uma unica raiz induzida pela relacao pai-filho da arvore materializada.

## Consequencia para o estudo multi-anchor

Este teorema explica os surveys em que havia muitos eventos multi-owner, mas
nenhum evento com varias ancoras no workspace geral. O fenomeno correto e:

```text
multi-owner em C
  nao implica
multi-anchor em W_general^-(C,b).
```

Assim, a preocupacao principal deixa de ser uma recomposicao por floresta de
varias raizes locais. A proxima obrigacao teorica e mostrar que a mutacao
local pode usar essa raiz unica como fechamento estrutural, mesmo quando a
implementacao atual escolhe apenas um owner extremal como `nodeCa`.

Em termos de codigo, a ponte natural e:

```text
closureRoot(C,b) = unico elemento de A_general(C,b).
```

O passo seguinte deve comparar `closureRoot(C,b)` com o `nodeCa` escolhido
pelo contrato atual e decidir se `updateTreeWithExplicitRegion` deve continuar
ancorado no owner extremal ou se deve passar a fechar a recomposicao pelo
representante unico de `K_a^-(C)`.

## Diagnostico `closureRoot` versus `nodeCa`

O survey experimental agora mede explicitamente:

```text
closureRoot(C,b) = unico elemento de A_general(C,b)
nodeCa = owner extremal escolhido pelo contrato operacional
```

e separa os eventos em:

```text
closure_root_matches_nodeca_events
closure_root_differs_from_nodeca_events
multi_owner_closure_root_differs_from_nodeca_events
nodeca_reaches_closure_root_events
nodeca_to_closure_root_same_level_path_events
closure_root_differs_but_same_level_path_events
```

Nos tres padroes sinteticos 6x6, com regioes conectadas ate tamanho 6,
260 regioes por imagem e `b in {0,...,15}`:

```text
accepted: 8200
matched rebuild: 8200
closure_root_events: 8200
closure_root_matches_nodeca_events: 434
closure_root_differs_from_nodeca_events: 7766
multi_owner_closure_root_differs_from_nodeca_events: 6073
nodeca_reaches_closure_root_events: 8200
nodeca_to_closure_root_same_level_path_events: 8200
closure_root_differs_but_same_level_path_events: 7766
max_nodeca_to_closure_root_path_length: 1
```

Nas duas imagens pseudoaleatorias 8x8, sementes 23 e 911, regioes conectadas
ate tamanho 8, 340 regioes por imagem e `b in {0,...,19}`:

```text
accepted: 6780
matched rebuild: 6780
closure_root_events: 6780
closure_root_matches_nodeca_events: 419
closure_root_differs_from_nodeca_events: 6361
multi_owner_closure_root_differs_from_nodeca_events: 5134
nodeca_reaches_closure_root_events: 6780
nodeca_to_closure_root_same_level_path_events: 6780
closure_root_differs_but_same_level_path_events: 6361
max_nodeca_to_closure_root_path_length: 1
```

O primeiro caso encontrado pela busca automatica ja ocorre com uma regiao de
um unico pixel:

```text
FIRST_CLOSURE_ROOT_MISMATCH
tree=max
C={0}
b=10
a=9
nodeCa=34
closureRoot=9
owners=1
materialized=1
```

Portanto, `nodeCa == closureRoot` nao e uma hipotese necessaria para os
marcos atuais baterem rebuild. A interpretacao mais provavel e que `nodeCa`
funciona como ancora operacional inferior da cadeia local, enquanto
`closureRoot` e a raiz teorica superior do workspace geral. A proxima prova
deve explicar por que a recomposicao por buckets, seguida de absorcao de nos
vazios e canonicalizacao same-level, transporta corretamente a cadeia de
`nodeCa` ate `closureRoot`.

Essa ponte foi formalizada em:

```text
docs/problem1_nodeca_closure_root_absorption.md
```

Conclusao: nos marcos atuais, `nodeCa` sempre alcanca `closureRoot(C,b)` por
um caminho same-level, e todo caso `nodeCa != closureRoot(C,b)` e explicado
por essa contracao canonica.
