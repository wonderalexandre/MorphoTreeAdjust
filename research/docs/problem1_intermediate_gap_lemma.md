# Lema da lacuna intermediaria implicita no Problema 1

Este texto formaliza a conjectura que surgiu dos surveys do Problema 1:
lacunas intermediarias sem candidato local nao precisam ser materializadas como
nos canonicos da arvore atualizada.

O ponto principal e fechar a completude local do workspace no caso atualmente
implementado: apos materializar `C`, todo componente antigo de nivel `v` que
pode contribuir com proper part fora de `C` para o componente de limiar de `C`
e alcancado pela varredura de fronteira e, portanto, aparece em
`MergeNodes(v,-)`.

Este ainda nao e um certificado completo de promocao do updater experimental
para producao. Ele fecha a parte teorica da lacuna intermediaria sob o contrato
single-`nodeCa` usado pelo workspace atual; a prova final do updater completo
ainda deve incluir a recomposicao de arestas e a canonicalizacao same-level.

## Notacao

Considere uma imagem `I^- : P -> V` e uma regiao conexa `C subset P`.
No caso max-tree do Problema 1:

```text
a = min I^-(C)
b > max I^-(C)
I^+(x) = b, se x in C
I^+(x) = I^-(x), caso contrario
Gamma(a,b) = {a+1, ..., b}
```

No caso dual min-tree:

```text
a = max I^-(C)
b < min I^-(C)
I^+(x) = b, se x in C
I^+(x) = I^-(x), caso contrario
Gamma(a,b) = {a-1, ..., b}
```

Para uniformizar a escrita, defina o conjunto de limiar dual:

```text
X_v(I) = {x in P : I(x) >= v}, para max-tree
X_v(I) = {x in P : I(x) <= v}, para min-tree
```

Se `C subset X_v(I^+)`, seja:

```text
K_v^+(C) = componente conexo de X_v(I^+) que contem C.
```

Chamaremos `K_v^+(C)` de componente de limiar de `C` no nivel `v`.

Na notacao do paper, o objeto local `C_v^+` corresponde a esse componente de
limiar. A questao e se `C_v^+` aparece como no materializado da arvore
canonica `T^+` ou se permanece apenas como objeto de prova.

## Materializacao canonica

Dizemos que o nivel `v` e materializado em `T^+` quando existe um no `N_v` tal
que:

```text
eta(N_v) = v
support(N_v) = K_v^+(C)
```

Caso contrario, dizemos que `C_v^+` e implicito na arvore canonica.

Em uma arvore canonica, um no materializado deve representar uma mudanca real
de limiar. Equivalentemente, ele deve ter proper parts diretas no nivel do no,
ou deve separar filhos/proper parts de modo que seu suporte nao seja
redundante em relacao aos niveis vizinhos. Um no sem proper parts diretas e
sem mudanca canonica de suporte e removido pela canonicalizacao.

## Workspace local

Apos a materializacao local de `C`, o workspace experimental calcula, para
cada `v in Gamma(a,b)`:

```text
MergeNodes(v,-)
```

Esse conjunto e formado por nos antigos visitados a partir da fronteira de
`C`, subindo pelos ancestrais relevantes ate o intervalo de recomposicao.

Uma lacuna de candidato ocorre quando:

```text
MergeNodes(v,-) = vazio.
```

Ela e intermediaria quando:

```text
v != b.
```

## Hipoteses

O argumento usa as seguintes hipoteses e um invariante de workspace que sera
demonstrado em seguida.

H1. `C` e conexo e o evento satisfaz o Problema 1, ou seu dual min-tree.

H2. A materializacao local de `C` ja foi aplicada quando `C` corta proper
parts diretas de nos antigos. Depois desse passo, os pixels de `C` que
participam do update sao proper parts diretas de nos locais apropriados.

H3. A arvore antiga e canonica: todo no vivo tem proper parts diretas, e nao
ha arestas pai-filho redundantes de mesmo nivel.

H4. Completude do workspace no intervalo afetado: se um no antigo de altitude
`v` tem proper part fora de `C` e contribui para `K_v^+(C) \\ C`, entao esse
no aparece em `MergeNodes(v,-)`. Este item e provado abaixo a partir da
varredura de fronteira e da ancestralidade na arvore de componentes.

H5. A unica alteracao de nivel da imagem ocorre em `C`, e todo pixel de `C`
recebe o nivel `b`.

## Fechamento de H4

Nesta secao consideramos o caso max-tree. O caso min-tree e obtido por
dualidade, trocando `>=` por `<=` e invertendo a ordem dos niveis.

No paper, os nos afetados sao definidos primeiro por conectividade e por
existencia no intervalo:

```text
MergeNodes(C,-)       = {X in Nodes^- : X union C e conexo}
MergeNodes([a,b],-)   = {X in Nodes^- : intervalle(X) intersect [a,b] != vazio}
MergeNodes^-          = MergeNodes(C,-) intersect MergeNodes([a,b],-)
```

Ha entao duas leituras por nivel que precisam ser separadas:

```text
MergeNodes_threshold(v,-) = MergeNodes^- intersect Nodes_v^-
MergeNodes_valuation(v,-) = {X in MergeNodes^- : Valuation^-(X) = v}
```

O objeto `C_v^+` e definido pelo componente de limiar, isto e, pela leitura
`threshold`. Ja o workspace experimental `MergeNodes(v,-)` e o bucket por
valuation: ele registra cada no antigo somente no nivel de sua propria
`Valuation^-`.

O criterio de materializacao real usa a leitura por valuation:

```text
Intervalle(C) =
  {v in [a,b-1] : MergeNodes_proper(v,-) != vazio} union {b},

MergeNodes_proper(v,-) =
  {X in MergeNodes_valuation(v,-) : Proper^-(X) \ C != vazio}.
```

Portanto, para um nivel intermediario `v != b`, so existe no canonico novo em
`v` quando algum no antigo de valuation `v`, conectado a `C`, ainda possui
proper part fora de `C`.

O workspace experimental implementa a mesma cobertura de forma local:

```text
1. percorre a fronteira grafica de C;
2. toma propernode^-(q) para cada pixel q adjacente a C e fora de C;
3. descarta sementes fora do lado correto de a;
4. sobe por Hasse^- enquanto o nivel ainda pertence ao intervalo de C;
5. registra cada no visitado no bucket de sua propria valuation.
```

### Lema de fronteira

Seja `X in Nodes^-`, com `Valuation^-(X)=v` e `v in [a,b]`, tal que:

```text
Proper^-(X) \ C != vazio
X union C e conexo em X_v(I^+)
```

Entao existe um pixel `q notin C`, adjacente a algum pixel de `C`, que pertence
ao suporte de `X`.

Prova. Como `X union C` e conexo em `X_v(I^+)`, existe um caminho de limiar
entre `C` e um pixel de `X \ C`. Tome a ultima saida desse caminho a partir de
`C` antes de atingir `X \ C`. O primeiro pixel fora de `C` depois dessa saida
e adjacente a `C`. Fora de `C`, a imagem nao mudou, logo esse trecho do
caminho pertence a `X_v(I^-)`. Como `X` e o componente antigo de nivel `v`
atingido por esse trecho, esse pixel de fronteira pertence ao suporte de `X`.

### Lema de ancestralidade

Se `q` pertence ao suporte de um no antigo `X in Nodes^-`, entao
`propernode^-(q)` e descendente de `X` ou igual a `X`. Assim, a subida por
`Hasse^-` a partir de `propernode^-(q)` visita `X`.

Prova. Em uma arvore de componentes canonica, os nos que contem um pixel fixo
formam uma cadeia totalmente ordenada por inclusao. O menor componente que
contem `q` e `propernode^-(q)`. Qualquer outro componente antigo que contem
`q`, em particular `X`, aparece como ancestral nessa cadeia. Portanto, subir
por `Hasse^-` a partir de `propernode^-(q)` alcanca `X`.

### Corolario de completude do workspace

Sob H1--H3 e H5, para todo `v in [a,b]`, se:

```text
Valuation^-(X) = v
Proper^-(X) \ C != vazio
X union C e conexo em X_v(I^+)
```

entao `X` aparece no bucket experimental `MergeNodes(v,-)`.

Prova. Pelo lema de fronteira, existe `q notin C` adjacente a `C` no suporte
de `X`. O workspace examina esse pixel de fronteira. Como `q in X_v(I^-)`,
`propernode^-(q)` tem valuation pelo menos `v`; como `v >= a`, a semente passa
pelo filtro `seedAltitude >= a`. Pelo lema de ancestralidade, a subida por
`Hasse^-` visita `X` antes de sair do intervalo de `C`, pois
`Valuation^-(X)=v >= a`. Como `v <= b`, `X` e registrado no bucket
`MergeNodes(v,-)`. Se `X` ja tiver sido marcado por outra semente, ele ja foi
registrado antes; a marcacao evita duplicacao, nao a cobertura.

Este corolario e exatamente H4 para o workspace single-`nodeCa`.

## Lema

Lema. Sob H1--H5, seja `v in Gamma(a,b)` com `v != b`. Se:

```text
MergeNodes(v,-) = vazio,
```

entao `C_v^+ = K_v^+(C)` nao e materializado como no canonico de `T^+`.
Equivalentemente, `C_v^+` e implicito na arvore canonica atualizada.

## Prova do lema

1. Pelo corolario de completude do workspace, o bucket experimental
`MergeNodes(v,-)` contem todo no antigo de valuation `v` que e conectado a
`C` e possui proper part fora de `C`.

2. Se `MergeNodes(v,-) = vazio`, entao necessariamente:

```text
MergeNodes_proper(v,-) = vazio.
```

3. Pelo criterio do paper:

```text
Intervalle(C) =
  {u in [a,b-1] : MergeNodes_proper(u,-) != vazio} union {b}.
```

Como `v != b` e `MergeNodes_proper(v,-) = vazio`, temos `v notin Intervalle(C)`.

4. A propriedade do paper que reduz os nos recompostos estabelece que a arvore
canonica precisa materializar apenas:

```text
{MergedNode(u,+) : u in Intervalle(C)}.
```

Logo, `MergedNode(v,+)` nao aparece como no canonico novo.

5. O componente de limiar `C_v^+ = K_v^+(C)` continua existindo como objeto de
prova: ele descreve o componente conexo de `X_v(I^+)` que contem `C`. Mas, como
`v notin Intervalle(C)`, esse componente coincide com o componente de algum
nivel vizinho materializado. Materializa-lo em `v` criaria uma cadeia
same-level ou um no sem proper part direta fora de `C`, que a canonicalizacao
remove.

6. Portanto, `C_v^+` e implicito na arvore canonica atualizada.

## Hasse+ com niveis implicitos

A recomposicao de arestas pode ser escrita diretamente sobre o conjunto
esparso de niveis materializados:

```text
S = Intervalle(C).
```

No caso max-tree, para qualquer `t in [a,b]`, defina:

```text
rep(t) = min {u in S : u >= t}.
```

Esse representante existe porque `b in S`. Se `t in S`, entao
`rep(t)=t`. Se `t notin S`, entao `C_t^+ = C_rep(t)^+`; isto e, o componente
de limiar em `t` existe como objeto de prova, mas sua representacao canonica e
o no materializado no proximo nivel de `S`.

Tambem defina, para `u in S`:

```text
pred_S(u) = max {w in S : w < u}, se esse conjunto nao e vazio;
pred_S(u) = alpha, caso contrario,
```

onde `alpha` e o primeiro nivel estavel abaixo de `a`, como na proposicao de
arestas do paper.

Com essa notacao, as arestas entre os novos nos recompostos ficam:

```text
MergedNode(u,+) Hasse+ MergedNode(pred_S(u),+), se pred_S(u) != alpha
MergedNode(u,+) Hasse+ MergedNode(alpha,+),     caso contrario.
```

As arestas que saem de nos antigos preservados tambem usam `rep`. Seja
`X Hasse- Y`, com `X` preservado em `T+`.

```text
1. Se X e uma fronteira acima de b:
      X Hasse+ MergedNode(b,+).

2. Se Y pertence a MergeNodes_valuation(t,-), com t in [a,b]:
      X Hasse+ MergedNode(rep(t),+).

3. Se Y nao pertence ao intervalo afetado:
      X Hasse+ Y.
```

Essa e exatamente a forma esparsa das proposicoes de arestas do paper. A forma
densa diz que `X` deve apontar para `MergedNode(t,+)`. Quando `t notin S`, a
propriedade `C_t^+ = C_rep(t)^+` troca esse alvo pelo representante canonico
sem alterar suporte, proper parts ou reconstrucao.

No caso min-tree dual, a mesma regra vale trocando a ordem:

```text
rep(t) = max {u in S : u <= t}
pred_S(u) = min {w in S : w > u}, se existir; caso contrario alpha.
```

Assim, niveis implicitos nao exigem arestas proprias em `Hasse+`: toda aresta
que apontaria para um `C_t^+` implicito aponta diretamente para seu
representante materializado.

## Canonicalizacao same-level

A implementacao experimental pode construir uma arvore intermediaria
`T~+` que ja reconstrui `I+`, mas ainda nao esta canonica. Isso acontece
quando a recomposicao local reaproveita candidatos existentes ou materializa
temporariamente um nivel que depois se mostra redundante. O sintoma estrutural
e uma aresta pai-filho de mesmo nivel:

```text
X Hasse~+ Y
eta(X) = eta(Y).
```

Defina a contracao:

```text
Proper(X)   := Proper(X) union Proper(Y)
Children(X) := Children(X) union Children(Y)
remove Y
```

Proposicao. Contrair uma aresta pai-filho de mesmo nivel preserva todos os
conjuntos de limiar da imagem reconstruida.

Prova. Como `eta(X)=eta(Y)`, os dois nos entram e saem dos mesmos conjuntos
`X_v(I)` para qualquer limiar `v`. Substituir a cadeia `X -> Y` pelo suporte
uniao nao cria nem remove componente conexo em nenhum limiar: os pixels de
`Y`, os filhos de `Y` e os proper parts de `Y` continuam no mesmo componente
representado por `X`. Portanto, a imagem reconstruida e os componentes
canonicos de limiar sao preservados.

Corolario. A contracao iterada de todas as arestas pai-filho de mesmo nivel
produz o quociente canonico de `T~+` em relacao aos niveis. O processo termina
porque cada contracao remove um no vivo. O resultado e independente da ordem
das contracoes no sentido estrutural: cada cadeia maximal conectada por
arestas same-level e substituida por um unico no com o suporte agregado.

Logo, se a recomposicao esparsa produz `T~+` com os mesmos componentes de
limiar de `I+`, entao:

```text
canonicalize_same_level(T~+) = T+ canonica.
```

Essa proposicao explica por que a etapa experimental
`collapseSameLevelParentChildEdges` nao muda a semantica do update: ela apenas
remove materializacoes locais que nao pertencem a `Intervalle(C)` ou arestas
redundantes criadas pela reutilizacao de candidatos.

## Corolario algoritmico

Sob as mesmas hipoteses, o algoritmo local nao precisa criar um no
intermediario para todo `v in Gamma(a,b)`.

O caminho suficiente e:

```text
1. materializar C quando houver owners parciais;
2. construir MergeNodes(v,-) para os niveis afetados;
3. criar explicitamente o no alvo b quando MergeNodes(b,-) estiver vazio;
4. recompor apenas os buckets nao vazios;
5. canonicalizar arestas pai-filho de mesmo nivel.
```

O unico caso que exigiria criar um no intermediario real seria um
contraexemplo do tipo:

```text
v != b
MergeNodes(v,-) = vazio
K_v^+(C) aparece como no canonico no rebuild
```

O diagnostico correspondente procura exatamente esse contraexemplo.

## Evidencia experimental

A busca automatica nao encontrou contraexemplo nos marcos atuais:

```text
sintetico 8x8:
  accepted=6780
  missing_intermediate_level_candidates=8190
  NO_MATERIALIZED_INTERMEDIATE_GAP

ICDAR 64x64 controlado:
  accepted=1017
  missing_intermediate_level_candidates=110359
  NO_MATERIALIZED_INTERMEDIATE_GAP

ICDAR 100x2 region16:
  accepted=762954
  missing_intermediate_level_candidates=57940973
  NO_MATERIALIZED_INTERMEDIATE_GAP
```

Esses resultados nao substituem a prova, mas validam a leitura operacional do
lema: nos eventos testados, toda lacuna intermediaria sem candidato se
comportou como objeto de prova, nao como no obrigatorio da arvore canonica.

## Proximas obrigacoes teoricas

O lema acima fecha a lacuna intermediaria para o workspace single-`nodeCa`,
apos materializacao local de `C`. As secoes seguintes tambem fixam a forma
esparsa de `Hasse+` e a canonicalizacao same-level como quociente canonico.
Para transformar isso em uma prova completa do updater experimental alem do
escopo single-`nodeCa`, ainda falta fechar a extensao de regiao geral
registrada em `research/docs/problem1_general_region_multi_anchor.md`:

```text
1. a extensao alem do contrato single-nodeCa, caso que ainda pertence ao plano
   de regiao geral/multiplos anchors.
```

As obrigacoes internas ao escopo single-`nodeCa` foram movidas para
`research/docs/problem1_explicit_region_mutation_proof.md`: equivalencia entre
`S_code` e `I_C^star`, recomposicao mutavel, quocientes por nos vazios,
canonicalizacao same-level e forma dual min-tree.

Depois dessa extensao, a recomposicao esparsa por buckets nao vazios fica
alinhada com a arvore canonica de `I^+` como teorema completo do updater local
para regioes gerais.
