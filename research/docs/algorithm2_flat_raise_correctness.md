# Corretude local do Algoritmo 2 restrito

Este texto formaliza o subcaso implementado por
`Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner`.

O resultado cobre apenas o contrato operacional descrito em
[`algorithm2_flat_raise_contract.md`](./algorithm2_flat_raise_contract.md). Ele
nao prova o problema geral de atualizacao local.

## Notacao

Seja `f : E -> V` uma imagem em um dominio discreto finito `E`, com uma
adjacencia fixa. Seja `T_f` a max-tree canonica de `f`. Para um no `Q`, denote
por `rho(Q)` o conjunto de proper parts diretas de `Q`, e por `level(Q)` sua
altitude.

Considere uma regiao `P subset E` e dois niveis `a,b in V`, com `b > a`.
Defina a imagem atualizada `f'` por

```text
f'(x) = b     se x in P,
f'(x) = f(x) caso contrario.
```

O algoritmo trabalha localmente em `T_f` e produz uma arvore `T'`.

## Teorema

Sob as hipoteses abaixo, `T'` e equivalente a max-tree canonica `T_{f'}` que
seria obtida por rebuild global de `f'`.

Hipoteses:

- `(H1)` `T_f` e uma max-tree.
- `(H2)` `P` e nao vazio.
- `(H3)` todo pixel de `P` pertence ao dominio `E` e aparece uma unica vez em
  `P`.
- `(H4)` `P` e conexo para a adjacencia usada por `T_f`.
- `(H5)` para todo `x in P`, `f(x)=a`.
- `(H6)` `b > a`.
- `(H7)` existe um unico no vivo `Q_a` tal que todo `x in P` tem `Q_a` como
  menor componente que o contem.
- `(H8)` `level(Q_a)=a`.
- `(H9)` todo `x in P` pertence a `rho(Q_a)`, isto e, e proper part direta de
  `Q_a`.
- `(H10)` `P = rho(Q_a)` ou `P` e uma parte propria de `rho(Q_a)`.

Conclusao:

```text
reconstruct(T') = f'
```

e a estrutura local de `T'`, apos contracao de arestas de mesma altitude, e a
mesma estrutura canonica que o rebuild global criaria dentro da regiao afetada.

## Hipoteses, guardas e testes

| Hipotese | Guarda no codigo | Evidencia de teste |
| --- | --- | --- |
| `(H1)` max-tree | `checkAlgorithm2FlatRaiseEligibility`: `maxTree.isMaxtree()` | `testAlgorithm2FlatRaiseRejectsMinTreeBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H2)` `P` nao vazio | `change.pixelsC.empty()` | `testAlgorithm2FlatRaiseRejectsEmptyRegionBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H3)` pixels validos e sem duplicata | bounds check e vetor `seen` em `checkAlgorithm2FlatRaiseEligibility` | `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H4)` `P` conexo | `isConnectedRegion` quando `assumeConnected=false` | `testAlgorithm2FlatRaiseRejectsDisconnectedRegionBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H5)` nivel atual igual a `a` | altitude do owner de cada pixel | `testAlgorithm2FlatRaiseRejectsLevelMismatchBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H6)` `b > a` | `change.b <= change.a` | `testAlgorithm2FlatRaiseRejectsNonRaiseBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H7)` unico `Q_a` vivo | owner comum dos pixels de `P` | `testAlgorithm2FlatRaiseRejectsMultipleSourceNodesBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `(H8)` `level(Q_a)=a` | altitude de `sourceNodeQa` | coberto por rejeicao de nivel e fixtures positivas |
| `(H9)` `P subset rho(Q_a)` | `isDirectProperPartOfNode` | fixtures positivas multi-pixel e regressao ICDAR |
| `(H10)` release ou preserve de `Q_a` | `eligibility.releasesQa` e comparacao `|P| == |rho(Q_a)|` | `testAlgorithm2FlatRaiseEligibilityAcceptsReleaseQa`, `testAlgorithm2FlatRaiseEligibilityAcceptsPreserveQa`, regressao ICDAR `release`/`preserve` |

As pos-condicoes estruturais sao verificadas por
`requireAlgorithm2FlatRaiseResultConsistent`, `requireNoEmptyLiveNodes` e
`requireStrictMaxTreeAltitudeOrder`.

## Objetos construidos

Como `P` muda de `a` para `b`, os unicos limiares que podem mudar sao:

```text
H^* = {h : a < h <= b}.
```

O algoritmo coleta dois tipos de fronteira local.

`Xi^\bullet_h` e o conjunto de nos de altitude `h`, com `a < h <= b`, que
tocam `P` pela adjacencia e devem ser absorvidos pela nova componente de nivel
`h`. No resultado do codigo, esses conjuntos aparecem em `xiBulletByLevel`.

`Xi^\circ` e o conjunto de nos de altitude estritamente maior que `b` que tocam
`P` e continuam sendo componentes mais altas. No resultado do codigo, esse
conjunto aparece em `xiCircleNodes`.

O algoritmo cria uma familia de niveis

```text
H = {b} union {h : Xi^\bullet_h nao e vazio}.
```

Para cada `h in H`, cria um novo no `Q_h`. O no `Q_b` recebe todos os pixels de
`P`. Os demais `Q_h` recompoem os limiares afetados abaixo de `b`.

## Lema 1: localidade

Todo componente de `T_f` que nao toca `P`, nem esta no caminho de ancestrais
locais de um componente que toca `P` dentro de `(a,b]`, mantem a mesma relacao
de inclusao em `T_{f'}`.

Prova. Fora de `P`, `f` e `f'` coincidem. Portanto, a unica forma de mudar um
componente de limiar superior e por meio da entrada de pixels de `P` nesse
limiar. Os componentes que podem ser conectados por essa entrada sao exatamente
os componentes adjacentes a `P` e seus ancestrais locais ate `b`. Esses sao os
componentes coletados pelo workspace como `Xi^\bullet` ou `Xi^\circ`. Logo, os
demais componentes permanecem inalterados.

## Lema 2: construcao correta de Q_b

No limiar `b`, os pixels de `P` passam a pertencer ao conjunto superior de
`f'`. Como `(H4)` garante que `P` e conexo, esses pixels formam uma unica nova
componente local. Componentes adjacentes de altitude maior que `b` continuam
sendo componentes estritamente mais altas, logo devem ficar abaixo de `Q_b` na
max-tree, e nao ser absorvidas por `Q_b`.

O algoritmo cria `Q_b`, move todos os pixels de `P` de `Q_a` para `Q_b` e
reanexa cada no de `Xi^\circ` sob `Q_b`. Assim, `Q_b` representa exatamente a
componente local que o rebuild global criaria no limiar `b`.

## Lema 3: recomposicao por nivel

Considere `h in H`, com `h < b`. Suponha que a cadeia criada para os niveis
maiores que `h` ja coincide com o rebuild global. No limiar `h`, a componente
local de `T_{f'}` e a uniao dessa cadeia mais alta com todos os componentes de
altitude `h` que passaram a se conectar a ela por meio de `P`. Esses
componentes sao precisamente `Xi^\bullet_h`.

O algoritmo cria `Q_h`, absorve todos os nos de `Xi^\bullet_h` nesse novo no e
anexa sob `Q_h` a componente criada no nivel imediatamente maior. Logo, por
inducao descendente nos niveis de `H`, todos os `Q_h` representam as
componentes locais que o rebuild global criaria.

## Lema 4: tratamento de Q_a

Ha dois casos.

Se `rho(Q_a)=P`, todos os proper parts diretos de `Q_a` foram movidos para
`Q_b`. A componente antiga de nivel `a` formada apenas por `P` deixa de ser
necessaria, pois `P` agora esta nos limiares superiores a `a`. O algoritmo
remove `Q_a` e coloca a raiz da recomposicao local em seu lugar.

Se `rho(Q_a) != P`, os proper parts de `Q_a` fora de `P` continuam no nivel
`a`. Portanto, a componente de nivel `a` deve continuar existindo. O algoritmo
preserva `Q_a`, remove apenas os pixels de `P` de `rho(Q_a)` e anexa a raiz da
recomposicao local abaixo de `Q_a`.

Assim, os dois ramos coincidem com a max-tree reconstruida de `f'`.

## Lema 5: canonicalizacao

A contracao final de arestas de mesma altitude nao altera a imagem
reconstruida. Ela apenas remove representacoes redundantes de uma mesma
componente no mesmo nivel. As validacoes finais garantem que:

- todo no vivo possui proper parts diretas;
- toda aresta da max-tree tem altitude estritamente crescente do pai para o
  filho;
- a reconstrucao de `T'` coincide com a imagem esperada nos testes.

Portanto, apos canonicalizacao, a estrutura local e canonica.

## Prova do teorema

Pelo Lema 1, todos os componentes fora da vizinhanca afetada permanecem
inalterados. Pelo Lema 2, `Q_b` representa corretamente a nova componente em
`b`. Pelo Lema 3, os demais `Q_h` recompoem todos os limiares afetados em
`(a,b)`. Pelo Lema 4, o antigo `Q_a` e removido ou preservado exatamente como
no rebuild global. Pelo Lema 5, a estrutura final fica canonica sem alterar a
imagem reconstruida.

Logo, sob `(H1)--(H10)`, `updateMaxTreeFlatRaiseSingleOwner` produz uma arvore
equivalente a `T_{f'}` no subcaso restrito.

## Limites

O teorema nao cobre:

- `P` distribuido em multiplos nos `Q_a`;
- `P` que exige materializacao para virar proper parts diretas de um unico
  `Q_a`;
- split em max-tree, isto e, `b < a`;
- versao dual para min-tree;
- manutencao incremental de atributos;
- integracao no CASF de producao.

Esses casos continuam como extensoes teoricas separadas.
