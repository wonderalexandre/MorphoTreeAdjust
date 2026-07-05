# Decisao de estabilizacao do Algoritmo 2 restrito

Data: 2026-06-04.

Branch: `theory/algorithm-alignment`.

Commit base da implementacao: `eb462e5 Stabilize restricted Algorithm 2 flat raise`.

## Decisao

Adotar a Opcao A do plano de promocao:

```text
manter em debug/opt-in, mas documentar como subcaso restrito estabilizado.
```

O algoritmo nao deve ser movido para uma API publica nem integrado ao CASF de
producao neste momento.

## Escopo estabilizado

O escopo estabilizado e exatamente:

- max-tree;
- atualizacao de aumento de nivel, `b > a`;
- regiao `P` plana, conexa, nao vazia, sem pixels duplicados e dentro do
  dominio;
- todos os pixels de `P` no nivel atual `a`;
- um unico no vivo `Q_a`;
- todos os pixels de `P` sao proper parts diretas de `Q_a`;
- `P = rho(Q_a)` ou `P` e parte propria de `rho(Q_a)`;
- atributo incremental validado no backend opt-in: `AREA`.

O gate operacional do escopo e
`checkAlgorithm2FlatRaiseEligibility` em
`research/include/morphoTreeAdjustResearch/debug/Algorithm2FlatRaiseEligibility.hpp`.

## Evidencia usada na decisao

Corretude local:

- contrato operacional em `algorithm2_flat_raise_contract.md`;
- prova local em `algorithm2_flat_raise_correctness.md`;
- matriz hipotese -> guarda -> teste;
- testes positivos, negativos e randomizados em
  `debug_region_merge_updater_unit_tests`.

Regressao automatizavel:

- `algorithm2_flat_raise_contract_regression` registrado no `ctest`;
- eventos reais `preserve Q_a` e `release Q_a` em imagem versionada;
- falha automatica com `--require-events`, `--fail-on-rejected` e
  `--fail-on-mismatch`.

Regressao ICDAR 1080p:

```text
case        images  events  rejected  reconstruction_mismatch  area_mismatch
preserve       100     500         0                        0              0
release        100     500         0                        0              0
```

CASF de producao:

- `dynamic_component_tree_casf_unit_tests`: passou;
- `dual_min_max_tree_incremental_filter_unit_tests`: passou;
- `public_core_api_unit_tests`: passou;
- `benchmark_component_tree_casf_smoke`: passou.

Resultados e comandos estao registrados em:

- `casf_performance_baseline.md`.

## Fora do escopo

Continuam fora desta decisao:

- `P` distribuido em multiplos nos `Q_a`;
- materializacao como parte do algoritmo provado;
- split em max-tree, `b < a`;
- versao dual para min-tree;
- atributos diferentes de `AREA`;
- filtros de atributo nao crescente;
- integracao no CASF de producao;
- substituicao de `updateTree` ou do caminho subtree-specialized atual.

## Justificativa

A evidencia atual e suficiente para classificar o metodo como backend restrito
estabilizado em `debug/opt-in`: o contrato esta explicito, os guardas rejeitam
entradas fora do escopo antes da mutacao, os dois ramos de `Q_a` foram
exercitados em testes e imagens reais, e o CASF de producao continuou passando.

Ainda nao ha justificativa para promover para API publica ou para o fluxo CASF
default, porque o metodo resolve apenas um subcaso local. Misturar essa
promocao com split, min-tree, materializacao geral ou atributos nao crescentes
abriria novamente o problema teorico geral.

## Condicoes para reabrir a promocao

Reavaliar a Opcao B somente se:

- a API experimental permanecer estavel depois de novos usos internos;
- os gates atuais continuarem passando;
- houver um consumidor opt-in claro fora dos testes;
- a manutencao de atributo do consumidor for `AREA` ou tiver prova e testes
  proprios.

Reavaliar integracao CASF somente se houver:

- equivalencia completa para a sequencia CASF alvo;
- fallback mensuravel e seguro;
- benchmark antes/depois sem regressao relevante no caminho default;
- decisao explicita separada desta nota.
