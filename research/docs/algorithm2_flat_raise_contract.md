# Algoritmo 2: contrato implementado para flat raise em max-tree

Este documento registra o escopo da implementacao experimental de
`Algorithm2FlatRegionMergeUpdater`. A implementacao e uma versao local e
construtiva do Algoritmo 2 do texto JMIV para o caso de aumento de nivel em
uma max-tree.

A prova local curta deste contrato esta em
[`algorithm2_flat_raise_correctness.md`](./algorithm2_flat_raise_correctness.md).

## Status de estabilizacao

Este contrato esta estabilizado como subcaso restrito em `debug/opt-in`.
Isso significa:

- ha API explicita de elegibilidade antes da mutacao da arvore;
- ha prova local para o caso restrito;
- ha testes unitarios positivos e negativos;
- ha regressao automatizavel via `ctest`;
- ha evidencia ICDAR 1080p para os ramos `rho(Q_a)=P` e `rho(Q_a)!=P`;
- o caminho CASF de producao permanece inalterado.

Esta classificacao nao torna o metodo uma solucao geral de atualizacao local.
Ela apenas estabiliza o subcaso descrito abaixo.

## Escopo de atributo

`Algorithm2FlatRegionMergeUpdater` e uma rotina estrutural: ela atualiza a
max-tree e retorna diagnosticos da recomposicao local. A manutencao de atributo
fica em uma camada separada.

No backend debug/opt-in `RegionMergeUpdater`, o uso do Algoritmo 2 em fluxo de
atualizacao fica restrito ao atributo `AREA`. Atributos `bbox_*` e atributos
radiometricos ficam fora do MVP e devem ser rejeitados nesse backend, ou
tratados por outro caminho sem usar a recomposicao local do Algoritmo 2.

## Entrada

A entrada e uma mudanca local de nivel

```text
P, a, b, com b > a.
```

No codigo, ela e representada por `RegionLevelChange<AltitudeType>`:

```cpp
RegionLevelChange<AltitudeType> change;
change.pixelsC = P;
change.a = a;
change.b = b;
change.direction = RegionUpdateDirection::Raise;
```

Apesar do nome historico `pixelsC`, neste contrato a regiao corresponde ao
conjunto `P` do Algoritmo 2.

## Hipoteses aceitas

A implementacao atual assume:

- a arvore atual e uma max-tree;
- `P` e nao vazio;
- os pixels de `P` pertencem ao dominio e nao aparecem duplicados;
- `P` e conexo para a adjacencia fornecida;
- todo pixel de `P` tem nivel atual `a`;
- `b > a`;
- existe um unico no `Q_a` tal que todo pixel de `P` e proper part direta de
  `Q_a`;
- `Q_a` tem altitude `a`;
- `P` pode ser todo `rho(Q_a)` ou apenas uma parte de `rho(Q_a)`.

A API estrita `updateMaxTreeFlatRaiseSingleOwner` verifica essas hipoteses
antes de executar a recomposicao. Se alguma delas falhar, a chamada deve
rejeitar a entrada sem alterar a arvore. O nome antigo
`updateMaxTreeFlatRaise` permanece apenas como wrapper de compatibilidade para
esse mesmo contrato restrito.

A implementacao cobre explicitamente os dois ramos finais do Algoritmo 2:

```text
rho(Q_a) = P      -> Q_a e removido;
rho(Q_a) != P     -> Q_a e preservado e P e removido de rho(Q_a).
```

## API de elegibilidade

O gate reutilizavel do contrato e:

```cpp
debug::problem2::checkAlgorithm2FlatRaiseEligibility(maxTree, adjacency, change);
debug::problem2::requireAlgorithm2FlatRaiseEligibility(maxTree, adjacency, change);
```

Arquivo:

```text
research/include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseEligibility.hpp
```

O resultado `Algorithm2FlatRaiseEligibility` informa se o evento foi aceito,
qual foi a razao de rejeicao, qual e `Q_a`, os niveis `a,b`, o tamanho de `P`,
o numero de proper parts diretas de `Q_a`, o numero de filhos de `Q_a`, o
numero de vizinhos de fronteira acima de `a`, e se o evento libera `Q_a`.

## Matriz de rastreabilidade

| Hipotese/decisao | Guarda principal | Evidencia |
| --- | --- | --- |
| max-tree | `checkAlgorithm2FlatRaiseEligibility`: `maxTree.isMaxtree()` | `testAlgorithm2FlatRaiseRejectsMinTreeBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `P` nao vazio | `change.pixelsC.empty()` | `testAlgorithm2FlatRaiseRejectsEmptyRegionBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| pixels no dominio | bounds check em `change.pixelsC` | `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| pixels sem duplicata | vetor `seen` em `checkAlgorithm2FlatRaiseEligibility` | `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `b > a` | `change.b <= change.a` | `testAlgorithm2FlatRaiseRejectsNonRaiseBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `P` conexo | `isConnectedRegion` quando `assumeConnected=false` | `testAlgorithm2FlatRaiseRejectsDisconnectedRegionBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| todo pixel de `P` no nivel `a` | altitude do owner de cada pixel | `testAlgorithm2FlatRaiseRejectsLevelMismatchBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| unico `Q_a` | todos os pixels apontam para o mesmo owner vivo | `testAlgorithm2FlatRaiseRejectsMultipleSourceNodesBeforeMutation`, `testAlgorithm2FlatRaiseEligibilityRejectsInvalidContracts` |
| `level(Q_a)=a` | altitude de `sourceNodeQa` | fixtures positivas e rejeicao por nivel |
| `P subset rho(Q_a)` como proper parts diretas | `isDirectProperPartOfNode` para todo pixel de `P` | fixtures positivas multi-pixel e regressao ICDAR |
| `rho(Q_a)=P` | `eligibility.releasesQa` | `testAlgorithm2FlatRaiseEligibilityAcceptsReleaseQa`, regressao ICDAR `release` |
| `rho(Q_a)!=P` | `!eligibility.releasesQa` | `testAlgorithm2FlatRaiseEligibilityAcceptsPreserveQa`, regressao ICDAR `preserve` |
| atributo do MVP | `RegionMergeUpdater::updatePairWithAlgorithm2AreaOptIn` aceita apenas `AREA` | estrategia em `algorithm2_flat_raise_attribute_strategy.md` e testes do backend opt-in |

## Objetos do paper expostos pelo resultado

O resultado `Algorithm2FlatRaiseApplyResult` expõe os objetos principais do
Algoritmo 2:

| Notacao do paper | Campo no codigo |
| --- | --- |
| `Q_a` | `sourceNodeQa` |
| `Q_b` | `nodeQb` |
| niveis afetados | `affectedLevels` |
| `Xi^\bullet` | `xiBulletNodes`, `xiBulletByLevel` |
| `Xi^\circ` | `xiCircleNodes` |
| novos `Q_h` | `createdNodesByLevel` |
| raiz local recomposta | `recompositionRoot` |
| nos absorvidos | `absorbedNodes` |
| nos liberados | `releasedNodes` |

Tambem ha flags diagnosticas:

- `rhoQaEqualsPBeforeUpdate`;
- `releasedSourceNodeQa`;
- `preservedSourceNodeQa`;
- `materializedBeforeUpdate`.

O caminho `updateMaxTreeFlatRaiseAfterMaterializingPartialOwners` preserva o
experimento anterior de tentar materializar proper parts antes da recomposicao.
Ele nao e o contrato estrito do Algoritmo 2; serve apenas como infraestrutura
experimental para estudar casos que ainda nao estao na prova local.

## Comportamento local

A recomposicao executa os seguintes passos:

1. Valida `P`, `a`, `b` e localiza `Q_a`.
2. Constroi o workspace de merge com `RegionMergeWorkspace`.
3. Classifica os nos afetados por nivel em `Xi^\bullet`.
4. Classifica as fronteiras preservadas acima de `b` em `Xi^\circ`.
5. Cria `Q_b` e os novos `Q_h` para os niveis afetados.
6. Absorve os nos de `Xi^\bullet` nos novos `Q_h`.
7. Move os pixels de `P` de `Q_a` para `Q_b`.
8. Reanexa `Xi^\circ` sob `Q_b`.
9. Encadeia os novos nos por nivel, em ordem compativel com a max-tree.
10. Remove ou preserva `Q_a` conforme `rho(Q_a)=P` ou `rho(Q_a)!=P`.
11. Canonicaliza arestas de mesma altitude e valida a estrutura final.

## Pos-condicoes verificadas

Os testes atuais verificam:

- a imagem reconstruida e igual ao rebuild global;
- todos os nos vivos possuem proper parts diretas;
- toda aresta da max-tree respeita altitude estritamente crescente do pai para
  o filho;
- `Q_b` e criado no nivel `b`;
- `Xi^\bullet` e absorvido nos novos nos do mesmo nivel;
- `Xi^\circ` e reanexado sob `Q_b`;
- `Q_a` e removido quando `rho(Q_a)=P`;
- `Q_a` e preservado quando `rho(Q_a)!=P`;
- entradas fora do contrato sao rejeitadas antes de alterar a reconstrucao da
  arvore;
- `P` com multiplos pixels conexos e testado explicitamente;
- fixtures 2D cobrem `Q_a` nao raiz, multiplas fronteiras `Xi^\circ` e
  multiplos niveis intermediarios em `Xi^\bullet`;
- uma fixture 2D cobre `Xi^\bullet` e `Xi^\circ` simultaneamente;
- um teste randomizado deterministico compara varias instancias pequenas contra
  rebuild.
- outro teste randomizado deterministico alterna os ramos
  `rho(Q_a)=P` e `rho(Q_a)!=P`.
- um teste randomizado deterministico cobre grades `7x7` e `9x9` dentro do
  contrato restrito.

## Fora do contrato atual

Ainda nao estao cobertos como algoritmo completo:

- regiao `P` distribuida em multiplos nos `Q_a`;
- regiao que nao pode ser representada como proper parts diretas de um unico
  `Q_a` sem materializacao;
- caso de max-tree com `b < a`, isto e split;
- versao simetrica em min-tree;
- integracao no caminho CASF de producao;
- manutencao incremental de atributos associada a esta recomposicao.

Esses casos devem continuar como pesquisa/experimento separado ate haver
prova local e validacao contra rebuild.
