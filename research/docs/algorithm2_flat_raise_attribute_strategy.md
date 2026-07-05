# Estrategia de atributos para o Algoritmo 2 restrito

Este documento define a estrategia de atributos para uma futura integracao do
subcaso restrito do Algoritmo 2:

```cpp
Algorithm2FlatRegionMergeUpdater::updateMaxTreeFlatRaiseSingleOwner
```

A decisao principal e manter a recomposicao estrutural desacoplada dos
atributos. O algoritmo estrutural deve produzir uma arvore correta e um
resultado diagnostico. A manutencao dos atributos deve ser uma camada separada.

## Estado atual

O CASF de producao usa `DynamicAttributeComputer` e buffers externos por arvore.
O fluxo atual em `DualMinMaxTreeIncrementalFilter` possui:

- `preProcessing(nodeId, buffer)`;
- `mergeProcessing(parentId, childId, buffer)`;
- `postProcessing(nodeId, buffer)`;
- hooks opcionais para `onMoveProperParts`, `onMoveProperPart` e
  `onNodeRemoved`;
- marcacao de nos afetados por `markAttributeUpdate`;
- recomputacao local por `computeAttributeOnTreeNode`.

Os atributos atualmente expostos no CASF publicado sao:

| Atributo | Implementacao atual | Tipo de regra |
| --- | --- | --- |
| `AREA` | `DynamicAreaComputer` | aditivo simples |
| `BOX_WIDTH` | `DynamicBoundingBoxComputer` | resumo geometrico |
| `BOX_HEIGHT` | `DynamicBoundingBoxComputer` | resumo geometrico |
| `DIAGONAL_LENGTH` | `DynamicBoundingBoxComputer` | resumo geometrico |

Variancia, media e soma de nivel de cinza ainda nao fazem parte da API publica
de atributos deste reposititorio.

## Problema

O Algoritmo 2 restrito altera a max-tree por:

- criacao de `Q_b`;
- criacao de nos intermediarios `Q_h`;
- movimento dos pixels de `P` de `Q_a` para `Q_b`;
- absorcao dos nos em `Xi^\bullet`;
- reanexo dos nos em `Xi^\circ`;
- remocao ou preservacao de `Q_a`;
- contracao final de arestas de mesma altitude.

Essas operacoes afetam atributos por dois mecanismos:

1. mudanca de proper parts diretas;
2. mudanca de relacoes pai-filho.

Logo, nao basta atualizar apenas `Q_b`. Tambem podem mudar os valores de seus
ancestrais, do antigo pai de uma fronteira reanexada, do antigo pai de um no
absorvido e do antigo pai de `Q_a` quando `Q_a` e removido.

## Regra conservadora inicial

Antes de promover para CASF real, a primeira regra segura deve ser:

```text
estrutura local pelo Algoritmo 2 + recomputacao global do atributo da max-tree
```

Essa regra nao e a meta final de performance, mas separa corretude estrutural
de corretude de atributo. Ela permite usar o Algoritmo 2 como backend opt-in
sem confiar ainda em uma regra incremental incompleta de atributos.

Criterio de aceite dessa regra:

- a imagem reconstruida bate com rebuild;
- a max-tree atualizada localmente bate com o rebuild em atributos calculados
  globalmente;
- o CASF de producao continua sem usar esse caminho por default.

## Regra local proposta

Depois da regra conservadora, a manutencao local deve usar uma lista explicita
de nos sujos.

O resultado estrutural atual ja expoe:

- `sourceNodeQa`;
- `nodeQb`;
- `recompositionRoot`;
- `createdNodes`;
- `createdNodesByLevel`;
- `xiBulletNodes`;
- `xiCircleNodes`;
- `absorbedNodes`;
- `releasedNodes`;
- `releasedSourceNodeQa`;
- `preservedSourceNodeQa`.

Para uma recomputacao local completa de atributos, ainda falta registrar antes
da mutacao:

- pai antigo de cada no em `Xi^\circ`;
- pai antigo de cada no em `Xi^\bullet`;
- pai antigo de `Q_a`;
- pai antigo de cada no que sera contraido na canonicalizacao.

Com essa informacao, a camada de atributos pode construir:

```text
D0 = nos criados
   union {Q_b, recompositionRoot}
   union nos vivos reanexados
   union pais antigos dos nos reanexados
   union pais antigos dos nos absorvidos
   union Q_a, se preservado
   union pai antigo de Q_a, se Q_a foi removido
```

Depois disso, a regra local segura e:

1. descartar de `D0` os nos mortos;
2. expandir `D0` para todos os ancestrais vivos ate a raiz;
3. ordenar os nos sujos por profundidade, dos filhos para os pais;
4. executar `preProcessing`, `mergeProcessing` e `postProcessing` em cada no
   sujo;
5. comparar contra recomputacao global nos testes.

Essa regra e conservadora: pode recomputar mais nos do que o minimo, mas evita
perder antigos ancestrais que tiveram filhos removidos ou reanexados.

## AREA

`AREA` e o primeiro atributo recomendado para validacao.

Regra:

```text
area(Q) = |rho(Q)| + sum area(child(Q)).
```

Como `AREA` nao precisa de cache persistente, basta recomputar os nos sujos em
ordem bottom-up usando `DynamicAreaComputer`. A primeira versao dessa camada
esta implementada em
`research/include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseAreaAttributeUpdater.hpp`.
Ela recebe o resultado estrutural do Algoritmo 2, monta as sementes sujas
especificas desse algoritmo e delega a recomputacao para
`research/include/morphoTreeAdjustResearch/debug/AttributeDirtyRecomputer.hpp`.

Criterio de teste:

- rodar `updateMaxTreeFlatRaiseSingleOwner`;
- recomputar localmente `AREA` usando a dirty closure;
- recomputar `AREA` globalmente na mesma arvore atualizada;
- comparar todos os nos vivos da arvore atualizada.

## Bounding box fora do MVP

`BOX_WIDTH`, `BOX_HEIGHT` e `DIAGONAL_LENGTH` usam summaries persistentes. Eles
ja possuem hooks para movimento de proper parts, mas o Algoritmo 2 restrito
ainda nao chama esses hooks nem garante o crescimento/invalidacao correta dos
summaries internos quando novos nos sao criados ou liberados.

A decisao atual e deixar `bbox_*` fora do MVP do Algoritmo 2. Portanto:

- o backend experimental baseado no Algoritmo 2 deve rejeitar raise local para
  atributos diferentes de `AREA`;
- o CASF de producao com `bbox_*` continua usando os caminhos ja existentes;
- o suporte incremental local para `bbox_*` so deve ser reaberto em uma frente
  separada, se houver motivacao depois de consolidar `AREA`.

## Soma, media e variancia de nivel de cinza

Esses atributos ainda nao estao na API publica atual. Se forem adicionados, a
base recomendada e um summary aditivo:

```text
count(Q)  = numero de pixels no subtree de Q
sum(Q)    = soma dos niveis dos pixels no subtree de Q
sumSq(Q)  = soma dos quadrados dos niveis dos pixels no subtree de Q
mean(Q)   = sum(Q) / count(Q)
var(Q)    = sumSq(Q) / count(Q) - mean(Q)^2
```

Para uma max-tree apos o raise, os pixels de `P` passam a contribuir com nivel
`b`, nao `a`, nos nos que os contem. Portanto, qualquer atributo radiometrico
precisa usar os valores atualizados de `f'`, nao os valores antigos.

Regra recomendada:

- implementar um `DynamicGrayLevelMomentComputer` separado;
- manter summaries locais `{count, sum, sumSq}`;
- atualizar por dirty closure inicialmente;
- so depois otimizar com hooks finos para `onMoveProperPart`.

## Mudancas de codigo necessarias

Antes de ativar atributos locais no Algoritmo 2 para CASF, implementar:

1. Um resultado estrutural mais rico, registrando pais antigos de nos
   reanexados, absorvidos, removidos e contraidos. Feito para o caso atual do
   Algoritmo 2 restrito.
2. Um builder de dirty closure independente do atributo. Feito em
   `AttributeDirtyRecomputer.hpp`.
3. Um executor de recomputacao bottom-up usando `DynamicAttributeComputer`.
   Feito em `AttributeDirtyRecomputer.hpp`; validado inicialmente com
   `DynamicAreaComputer`. O mesmo helper tambem oferece recomputacao global
   bottom-up como fallback experimental quando nao houver regra local validada.
4. Testes para `AREA`. Feito com casos `Q_a` preservado e `Q_a` liberado.
5. Guarda explicita no backend experimental: Algoritmo 2 local apenas quando o
   atributo e `AREA`. Feito em `RegionMergeUpdater`.
6. Opcionalmente, novo computador de momentos de nivel de cinza em uma frente
   separada.

## Gate para CASF

O Algoritmo 2 restrito so deve ser usado em fluxo CASF depois que:

- `AREA` passar contra recomputacao global;
- `bbox_*` estiver explicitamente fora do backend Algoritmo 2;
- atributos radiometricos estiverem explicitamente fora do backend Algoritmo 2;
- benchmarks incluirem tempo estrutural e tempo de atributo separadamente.

## Decisao atual

A decisao atual e:

```text
MVP do Algoritmo 2 restrito = max-tree, flat raise `b > a`, contrato
single-owner/materialized, atributo `AREA`.
```

O estado atual tem dois caminhos para atributos:

- caminho local: dirty closure de ancestrais, validado para `AREA`;
- fallback global experimental: recomputacao bottom-up de todos os nos vivos
  usando `DynamicAttributeComputer`, disponivel como infraestrutura, mas nao
  usado para habilitar `bbox_*` no MVP.

O backend experimental `RegionMergeUpdater` deve aceitar o caminho Algoritmo 2
apenas quando `attribute == AREA`; atributos diferentes devem ser rejeitados ou
encaminhados a outro backend sem usar a recomposicao local do Algoritmo 2.

O uso recomendado do MVP esta descrito em
[`algorithm2_area_mvp_usage.md`](./algorithm2_area_mvp_usage.md).

O proximo passo pratico e consolidar benchmark e integracao opt-in para `AREA`,
sem abrir `bbox_*` nesta etapa.
