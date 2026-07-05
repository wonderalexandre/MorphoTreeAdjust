# Problema 1 para regioes gerais multi-anchor

Este texto inicia a extensao teorica do Problema 1 para regioes conexas `C`
que nao satisfazem o contrato single-`nodeCa` ja fechado em:

```text
research/docs/problem1_explicit_region_mutation_proof.md
```

O objetivo aqui nao e promover uma implementacao nova. O objetivo e separar o
que ja esta provado, enunciar o obstaculo multi-anchor e definir o programa de
estudo antes de alterar o caminho experimental.

## Ponto de partida

Considere uma imagem `I^- : P -> V`, uma arvore canonica `T^-` e uma regiao
conexa `C subset P`. O evento do Problema 1 max-tree e:

```text
a = min I^-(C)
b > max I^-(C)
I^+(x) = b, se x in C
I^+(x) = I^-(x), caso contrario
```

No dual min-tree:

```text
a = max I^-(C)
b < min I^-(C)
I^+(x) = b, se x in C
I^+(x) = I^-(x), caso contrario
```

Usaremos a notacao unificada:

```text
X_v(I) = {x in P : I(x) >= v}, para max-tree
X_v(I) = {x in P : I(x) <= v}, para min-tree

K_v^+(C) = componente conexo de X_v(I^+) que contem C.
```

O caso ja fechado assume, apos a materializacao local de `C`, que a regiao
afetada tem uma unica ancora estrutural `nodeCa`. Nesse caso, a prova mostra
que a recomposicao por buckets:

```text
S_code = I_C^star
Intervalle(C) subset S_code
b in S_code
```

seguida de absorcao de nos vazios e canonicalizacao same-level produz a
arvore canonica de `I^+`, no escopo single-`nodeCa`.

## Onde o caso geral difere

Depois da materializacao local de `C`, defina:

```text
O(C) = {propernode^-(x) : x in C}
```

isto e, os owners diretos dos pixels de `C` antes do update. O updater atual
escolhe um unico owner extremal:

```text
max-tree: owner de menor altitude em O(C)
min-tree: owner de maior altitude em O(C)
```

Esse owner funciona como `nodeCa` quando ele e suficiente para fechar
topologicamente todo o workspace afetado. O caso multi-anchor aparece quando a
regiao afetada possui mais de uma raiz local no intervalo ativo. Uma forma
operacional de registrar isso e definir o subgrafo afetado antigo:

```text
W^-(C,b) =
  nos antigos conectados a C e com altitude no intervalo ativo ate b,
  incluindo owners cortados pela materializacao e nos alcancados pela
  varredura de fronteira.
```

As ancoras locais sao:

```text
A(C,b) = {X in W^-(C,b) : parent^-(X) notin W^-(C,b)}
```

No contrato ja provado:

```text
|A(C,b)| = 1.
```

No caso geral:

```text
|A(C,b)| > 1.
```

Essa diferenca e estrutural. Nao basta escolher uma das ancoras por altitude,
pois a recomposicao final precisa recolocar a nova cadeia de componentes no
lugar de uma floresta local, nao apenas no lugar de uma subarvore antiga.

Observacao operacional importante: a primeira metrica implementada para
`A(C,b)` usava apenas o workspace esparso de `RegionMergeWorkspace`. Esse
workspace ainda era influenciado pela escolha de um unico `nodeCa`, portanto
ele era util como sinal de regressao, mas nao bastava como definicao geral.
A instrumentacao agora separa tres objetos:

```text
W_sparse^-(C,b):
  workspace usado pela recomposicao atual a partir de um nodeCa escolhido;

W_owner^-(C,b):
  ancestrais, dentro do intervalo ativo, de todos os owners de C;

W_general^-(C,b):
  W_owner^-(C,b) union W_sparse^-(C,b).
```

As ancoras relevantes para decidir se ha floresta local passam a ser:

```text
A_general(C,b) =
  {X in W_general^-(C,b) : parent^-(X) notin W_general^-(C,b)}.
```

Assim, `multiple_sparse_workspace_anchor_events` deve ser lido como
diagnostico historico do codigo atual, enquanto
`multiple_general_workspace_anchor_events` e o sinal mais proximo do problema
teorico.

## Problema geral a provar

Enunciado desejado, ainda como alvo:

```text
Entrada:
  T^- canonica;
  C conexo;
  evento max-tree ou min-tree do Problema 1;
  materializacao local de todos os owners cortados por C.

Construir localmente:
  T^+ canonica de I^+

usando apenas:
  buckets locais por valuation;
  fronteira grafica de C;
  criacao opcional do no alvo b;
  recomposicao local;
  absorcao de nos vazios;
  canonicalizacao same-level.
```

O novo ponto de prova e substituir o fechamento unico em `nodeCa` por uma
regra de fechamento em floresta:

```text
finalUnionNode + A(C,b) + fronteira externa do workspace
  -> mesma relacao Hasse+ obtida por rebuild.
```

## Grafo de incidencia das ancoras

Para entender a interacao entre ancoras, defina uma atribuicao local:

```text
owner_C(x) = propernode^-(x), para x in C.
```

Cada pixel de `C` pertence a algum owner materializado. Para cada ancora
`A_i in A(C,b)`, defina inicialmente:

```text
C_i = {x in C : owner_C(x) pertence ao ramo de A_i em W^-(C,b)}.
```

Os conjuntos `C_i` podem falhar em ser conexos. Esse e um dos primeiros
obstaculos teoricos: o Problema 1 assume `C` conexo, mas uma decomposicao por
ancora pode produzir partes desconexas.

Defina entao o grafo de incidencia:

```text
G_A(C) = (A(C,b), E_A)
```

com uma aresta `A_i -- A_j` quando existe uma adjacencia grafica
`p ~ q`, com:

```text
p in C_i
q in C_j
i != j.
```

Como `C` e conexo, `G_A(C)` deve capturar como as subregioes de ancoras
distintas sao conectadas pela propria regiao modificada. Quando `G_A(C)` tem
mais de um vertice, a criacao do componente alvo em `b` tende a unir ramos que
eram independentes em `T^-`.

## Estrategia 1: decomposicao sequencial

A primeira possibilidade a estudar e decompor:

```text
C = C_1 union ... union C_k
```

e aplicar o algoritmo single-`nodeCa` em cada parte. Para isso ser correto,
precisariamos de um teorema de comutacao:

```text
Update(C_i,b) seguido de Update(C_j,b)
  equivale a
Update(C_i union C_j,b)
```

para qualquer ordem admissivel.

Essa estrategia e atraente porque reutiliza a prova existente. Mas ela exige
hipoteses fortes:

```text
D1. cada C_i e conexo;
D2. cada C_i satisfaz o contrato single-nodeCa apos materializacao;
D3. o update de C_i nao altera os owners necessarios para validar C_j;
D4. as cadeias esparsas de S_code(C_i) e S_code(C_j) so se encontram no alvo b
    ou em uma regiao ja canonica comum;
D5. a ordem das atualizacoes nao muda a arvore canonica final.
```

Se alguma dessas condicoes falhar, a decomposicao sequencial deixa de ser uma
prova simples. Em particular, `D3` e suspeita: depois que uma parte de `C` vai
para `b`, ela pode virar conector para a parte seguinte e mudar o workspace
observado pelo segundo update.

Conjectura de decomposicao limitada:

```text
Se D1-D5 valem, entao a decomposicao sequencial e correta e o resultado
independe da ordem.
```

Essa conjectura deve ser testada primeiro por oracle, nao implementada como
caminho de producao.

## Estrategia 2: recomposicao simultanea

A segunda possibilidade e tratar `C` inteiro de uma vez. O workspace deixa de
ser uma cadeia fechada por um unico `nodeCa` e passa a ser uma floresta local:

```text
W^-(C,b), A(C,b), S_code(C,b).
```

O sweep por nivel continuaria global:

```text
max-tree: niveis de S_code em ordem decrescente a partir de b
min-tree: niveis de S_code em ordem crescente a partir de b
```

O alvo `b` continuaria unico, porque `C` e conexo e todos os pixels de `C`
recebem o mesmo nivel. Assim:

```text
K_b^+(C)
```

deve ser um unico componente alvo.

A dificuldade fica no fechamento inferior/superior do workspace. Em vez de:

```text
finalUnionNode -> nodeCa
```

precisamos recompor:

```text
finalUnionNode -> boundary(W^-(C,b))
```

onde `boundary(W^-(C,b))` contem os pais externos das ancoras removidas ou
preservadas. O teorema precisa dizer quais desses pais externos continuam
vivos, quais ancoras sao absorvidas e como filhos preservados devem ser
reanexados.

Conjectura simultanea:

```text
A recomposicao por um unico S_code global, seguida de fechamento por floresta
e canonicalizacao same-level, produz a mesma arvore canonica que o rebuild.
```

Essa estrategia parece mais proxima do enunciado geral do Problema 1, mas
exige uma nova regra mutavel de fechamento.

## Hipoteses simplificadoras candidatas

Antes de implementar, vale estudar hipoteses intermediarias. Elas podem gerar
teoremas uteis e benchmarks controlados.

H-MA1. `C_i` conexo para toda ancora.

```text
Permite testar decomposicao sequencial sem resolver subpartes desconexas.
```

H-MA2. `G_A(C)` e uma arvore.

```text
Evita ciclos de incidencia entre ancoras e torna a ordem de contracao mais
facil de especificar.
```

H-MA3. Todas as ancoras possuem a mesma altitude extremal.

```text
Reduz o fechamento multi-anchor a uma camada horizontal antes de conectar ao
alvo b.
```

H-MA4. Os workspaces das ancoras sao disjuntos fora de `C`.

```text
Evita que duas partes disputem o mesmo no antigo em buckets por valuation.
```

H-MA5. Existe um unico pai externo comum para `A(C,b)`.

```text
Substitui o fechamento em floresta por um fechamento unico em um ancestral
externo comum, possivelmente criado ou materializado.
```

Minha leitura inicial: H-MA1 e H-MA4 sao boas para testar decomposicao; H-MA5
e boa para aproximar a recomposicao simultanea sem resolver a floresta geral
de uma vez.

## Plano de oracle

Antes de mudar o updater, o oracle deve enumerar eventos multi-anchor pequenos
e responder perguntas estruturais.

A primeira instrumentacao desse plano esta em:

```text
research/include/morphoTreeAdjustResearch/debug/ExplicitRegionProblem1MultiAnchor.hpp
```

Ela foi usada com duas flags de diagnostico:

```text
--analyze-multi-anchor
--print-first-multi-anchor
--print-first-workspace-multi-anchor
```

Essa instrumentacao ainda nao prova a regra multi-anchor. Ela calcula sinais
estruturais para eventos aceitos pelo Problema 1:

```text
owners de C apos materializacao opcional;
owners extremais no nivel a;
componentes conexos de C por owner;
arestas de incidencia entre owners dentro de C;
raizes do workspace esparso usado pela recomposicao atual;
raizes do workspace por ancestrais de owners no intervalo ativo;
raizes do workspace geral W_general^-(C,b).
```

1. Classificacao de evento:

```text
calcular O(C)
materializar C
construir W^-(C,b)
calcular A(C,b)
se |A(C,b)| == 1: caso ja coberto
se |A(C,b)| > 1: caso multi-anchor
```

2. Diagnostico de decomposicao:

```text
calcular C_i por ancora
verificar conectividade de cada C_i
testar todas as ordens pequenas de Update(C_i,b)
comparar cada resultado com rebuild de Update(C,b)
registrar se a ordem altera a arvore final
```

3. Diagnostico de workspace simultaneo:

```text
construir S_code(C,b) global
comparar S_code(C,b) com union_i S_code(C_i,b)
detectar buckets compartilhados por mais de uma ancora
detectar pais externos distintos em boundary(W^-(C,b))
```

4. Busca de contraexemplos:

```text
primeiro caso em que decomposicao falha;
primeiro caso em que C_i e desconexo;
primeiro caso com buckets compartilhados;
primeiro caso com mais de um pai externo;
primeiro caso em que uma ordem sequencial bate rebuild e outra nao.
```

Esses contraexemplos devem virar figuras e lemas. Eles dirao se o caminho
correto e uma decomposicao com hipoteses ou uma recomposicao simultanea geral.

Exemplo inicial encontrado pela nova instrumentacao em um padrao sintetico
6x6:

```text
FIRST_MULTI_ANCHOR
tree=max
C={0,1,2}
b=6
a=0
owners=3
extremal_owners=2
owner_partition_components=3
owner_incidence_edges=2
sparse_workspace_anchors=1
materialized=1
```

Esse caso mostra que a regiao pode ter multiplos owners extremais mesmo
quando o workspace esparso ainda fecha em uma unica raiz. Portanto, o primeiro
subproblema experimental e distinguir:

```text
multi-owner em C
multi-anchor real no workspace
```

Eles nao sao automaticamente equivalentes.

Em seguida, a instrumentacao foi refinada para separar:

```text
multiple_owner_events:
  C possui mais de um owner apos materializacao opcional;

multi_owner_single_workspace_events:
  C possui mais de um owner, mas o workspace geral ainda tem uma unica raiz;

multiple_sparse_workspace_anchor_events:
  o workspace esparso possui mais de uma raiz local;

multiple_general_workspace_anchor_events:
  W_general^-(C,b) possui mais de uma raiz local.
```

Na varredura sintetica inicial com tres padroes 6x6, regioes conectadas
amostradas ate tamanho 6, 260 regioes por imagem e `b in {0,...,15}`:

```text
accepted: 8200
multiple_owner_events: 6295
multi_owner_single_workspace_events: 6295
multi_anchor_candidates: 581
multiple_extremal_owner_events: 573
disconnected_owner_partition_events: 8
multiple_sparse_workspace_anchor_events: 0
multiple_general_workspace_anchor_events: 0
max_sparse_workspace_anchors: 1
max_owner_interval_workspace_anchors: 1
max_general_workspace_anchors: 1
```

Portanto, nesse marco inicial, os sinais multi-owner sao frequentes, mas
nenhum evento exigiu uma floresta de workspace, nem pela metrica esparsa nem
pela metrica geral.

A busca foi entao ampliada para duas imagens pseudoaleatorias 8x8, sementes
fixas 23 e 911, regioes conectadas ate tamanho 8, 340 regioes por imagem e
`b in {0,...,19}`:

```text
accepted: 6780
multiple_owner_events: 5355
multi_owner_single_workspace_events: 5355
multi_anchor_candidates: 152
multiple_extremal_owner_events: 45
disconnected_owner_partition_events: 107
multiple_sparse_workspace_anchor_events: 0
multiple_general_workspace_anchor_events: 0
max_sparse_workspace_anchors: 1
max_owner_interval_workspace_anchors: 1
max_general_workspace_anchors: 1
NO_WORKSPACE_MULTI_ANCHOR
```

Esse resultado muda a prioridade teorica. A pergunta imediata deixa de ser
"como recompor uma floresta com varias ancoras encontradas pelo survey?" e
passa a ser:

```text
Conjectura de raiz unica local:
  para C conexo, apos materializacao local dos owners cortados por C,
  W_general^-(C,b) possui exatamente uma ancora local.
```

Se essa conjectura for verdadeira, os casos com multiplos owners e multiplos
owners extremais nao exigem necessariamente uma recomposicao por floresta. O
proximo teorema deve mostrar por que a conexidade de `C` e o fato de todos os
pixels de `C` irem para o mesmo nivel `b` forcariam a conexidade ancestral do
workspace geral. Se a conjectura falhar, o primeiro contraexemplo deve ser
gerado por `--print-first-workspace-multi-anchor` e virar o caso minimo para a
regra simultanea por floresta.

Esse teorema foi formalizado em:

```text
research/docs/problem1_general_workspace_unique_root.md
```

Conclusao parcial: sob as hipoteses do Problema 1, apos materializacao local
dos owners parciais de `C`, o workspace geral possui uma unica ancora:

```text
|A_general(C,b)| = 1.
```

A prova identifica essa ancora como o representante materializado de
`K_a^-(C)`, o componente conexo de `X_a(I^-)` que contem `C`. Portanto, os
casos com multiplos owners nao obrigam, por si so, uma recomposicao por
floresta. Eles ainda exigem uma prova de que a mutacao local pode fechar no
representante unico correto, e nao apenas em um owner extremal arbitrario.

O diagnostico `closureRoot` versus `nodeCa` mostrou que esses nos
frequentemente diferem, mas os marcos testados ainda batem rebuild:

```text
6x6 patterns:
  closure_root_differs_from_nodeca_events = 7766 / 8200
  matched rebuild = 8200 / 8200

8x8 PRNG:
  closure_root_differs_from_nodeca_events = 6361 / 6780
  matched rebuild = 6780 / 6780
```

Assim, a pergunta correta nao e se `nodeCa` coincide com a raiz teorica. A
pergunta e se a cadeia local iniciada em `nodeCa` sempre e absorvida e
canonicalizada ate a mesma estrutura que teria sido obtida fechando
diretamente em `closureRoot`.

Essa pergunta foi refinada em
`research/docs/problem1_nodeca_closure_root_absorption.md`. Nos marcos sinteticos
atuais, todo `nodeCa` alcanca `closureRoot` por um caminho same-level, com
comprimento maximo 1. Assim, a diferenca observada entre a ancora operacional
e a raiz teorica e explicada por uma contracao canonica local, nao por uma
nova floresta de fechamento.

O segundo diagnostico separa a completude dos buckets em uma condicao mais
fraca. Os buckets de fronteira nao contem todos
os nos do workspace owner-interval, mas nenhum caso testado deixou no ausente
sem explicacao. Todo no ausente caiu em uma destas classes:

```text
1. pertence ao caminho same-level nodeCa -> closureRoot(C,b);
2. possui apenas proper parts diretas contidas em C.
```

Portanto, a proxima prova nao deve tentar demonstrar a igualdade ingenua
entre workspace owner-interval e buckets de fronteira. Ela deve demonstrar que
todo no com uma proper part direta externa a `C`, fora do caminho de
fechamento, precisa aparecer em algum bucket `MergeNodes_valuation(v,-)`.
Essa formulacao foi registrada como lema em
`research/docs/problem1_frontier_bucket_completeness_lemma.md`.

Com esse lema, a raiz unica do workspace geral e a absorcao
`nodeCa -> closureRoot`, o criterio local suficiente para a correcao foi
integrado em `research/docs/problem1_integrated_local_correctness_theorem.md`.

## Criterio de avanco

O proximo marco teorico deve produzir uma das duas conclusoes:

```text
Conclusao A:
  A cadeia operacional nodeCa -> closureRoot e same-level e pode ser
  quotientada pela canonicalizacao final.

Conclusao B:
  Mesmo com fechamento resolvido, a mutacao atual ainda falha em casos gerais
  porque os buckets de merge ou filhos externos nao cobrem todos os ramos
  afetados.
```

Somente depois dessa decisao faz sentido implementar a generalizacao mutavel.
Enquanto isso, `updateTreeWithExplicitRegion` deve permanecer promovido apenas
no escopo ja provado e validado por rebuild.
