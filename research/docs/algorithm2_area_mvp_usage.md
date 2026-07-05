# Guia de uso do MVP AREA para o Algoritmo 2

Este guia descreve o uso recomendado do backend experimental opt-in baseado no
Algoritmo 2 restrito.

## Escopo

O MVP cobre somente:

- max-tree;
- flat raise, isto e, `b > a`;
- regiao `P` plana;
- um unico no fonte `Q_a`;
- todos os pixels de `P` sao proper parts diretas de `Q_a`;
- atributo `AREA`.

Ficam fora deste MVP:

- `bbox_*`;
- soma, media e variancia de nivel de cinza;
- split, isto e, `b < a`;
- multiplos nos fonte `Q_a`;
- integracao no CASF de producao.

## API opt-in

A entrada recomendada e:

```cpp
debug::region::RegionMergeAlgorithm2AreaOptInDiagnostics diagnostics;
std::vector<float> maxArea = /* AREA atual da max-tree */;

auto updated = debug::region::RegionMergeUpdater::updatePairWithAlgorithm2AreaOptIn(
    minTree,
    maxTree,
    adjacency,
    change,
    AREA,
    &maxArea,
    &diagnostics);
```

Comportamento:

- se `attribute == AREA` e o contrato do Algoritmo 2 for satisfeito, a max-tree
  e atualizada localmente, `maxArea` e atualizado por dirty closure e a min-tree
  e reconstruida a partir da imagem atualizada;
- se `attribute != AREA`, o backend nao tenta Algoritmo 2 e usa fallback por
  rebuild;
- se o contrato do Algoritmo 2 falhar, o backend usa fallback por rebuild.

O fallback reconstroi max-tree e min-tree a partir da imagem anterior mais a
mudanca `change`, preservando a corretude do resultado experimental.

## Diagnosticos

O struct `RegionMergeAlgorithm2AreaOptInDiagnostics` registra:

- `attemptedAlgorithm2`: houve tentativa de usar Algoritmo 2;
- `acceptedAlgorithm2`: a tentativa foi aceita;
- `rejectedByAttribute`: `attribute != AREA`;
- `rejectedByContract`: o contrato do Algoritmo 2 falhou;
- `fallbackUsed`: o resultado veio de rebuild;
- `areaAttributeUpdated`: o buffer `AREA` foi atualizado localmente;
- `createdNodes`, `absorbedNodes`, `releasedNodes`, `affectedLevels`;
- `xiBulletNodes`, `xiCircleNodes`;
- `areaDirtySeeds`, `areaDirtyClosure`, `areaRecomputedNodes`;
- `algorithm2StructuralMicros`, `areaAttributeMicros`, `primalRebuildMicros`,
  `fallbackRebuildMicros`.

## Regra pratica

Para experimentos do MVP:

```text
usar Algoritmo 2 somente quando diagnostics.acceptedAlgorithm2 == true
```

Para comparacoes de tempo, medir separadamente:

- tempo estrutural;
- tempo de `AREA` local;
- tempo de fallback/rebuild.

## Arquivos relevantes

- `research/include/morphoTreeAdjustResearch/debug/RegionMergeUpdater.hpp`
- `research/include/morphoTreeAdjustResearch/debug/Algorithm2FlatRegionMergeUpdater.hpp`
- `research/include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseAreaAttributeUpdater.hpp`
- `research/include/morphoTreeAdjustResearch/debug/AttributeDirtyRecomputer.hpp`
- `research/tests/debug_region_merge_updater_unit_tests.cpp`
