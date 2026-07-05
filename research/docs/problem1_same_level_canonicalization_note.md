# Nota: canonicalizacao same-level no Problema 1

Esta nota registra a interpretacao teorica da etapa experimental que colapsa
arestas pai-filho de mesmo nivel apos a atualizacao local do Problema 1.

## Contexto

No Problema 1, atualizamos uma arvore dual a partir de uma regiao conexa
explicita `C` e de um novo nivel `b`. Para a max-tree, a hipotese publicada e:

```text
C subset P conexo
a = min I-(C)
b > max I-(C)
I+(x) = b, se x em C
I+(x) = I-(x), caso contrario
```

A forma dual para a min-tree troca a direcao:

```text
b < min I-(C)
```

O oracle estrutural compara a arvore atualizada localmente com a arvore
reconstruida do zero a partir de `I+`. A comparacao nao depende dos ids dos
nos; ela usa altitude, suporte completo, proper parts diretos e relacao
pai-filho por suporte.

## Fenomeno observado

Apos materializacao local e recomposicao por varredura de niveis, a arvore
local pode conter uma aresta pai-filho com a mesma altitude. Exemplo observado:

```text
max-tree: parent_alt = 9, child_alt = 9
```

Essa arvore ainda reconstrui a imagem correta, mas nao e uma arvore canonica:
em uma max-tree canonica, uma aresta pai-filho deve ter altitude estritamente
crescente; em uma min-tree canonica, estritamente decrescente.

## Contracao canonica

Se um no `Y` e filho direto de `X` e ambos possuem a mesma altitude, entao a
aresta `X -> Y` nao representa uma mudanca de limiar. A contracao:

```text
Proper(X) := Proper(X) union Proper(Y)
Children(X) := Children(X) union Children(Y)
remove Y
```

preserva:

- a imagem reconstruida, pois `alt(X) = alt(Y)`;
- o suporte representado pela subarvore contraida;
- a particao dos pixels em proper parts diretos;
- a ordem canonica de niveis apos a remocao da aresta redundante.

Portanto, a contracao same-level deve ser interpretada como uma etapa de
canonicalizacao da recomposicao local, nao como uma mudanca semantica do
operador.

Mais precisamente, se a recomposicao local produz uma arvore intermediaria
`T~+`, a canonicalizacao same-level e o quociente de `T~+` pela relacao gerada
por arestas pai-filho de mesma altitude. Cada cadeia maximal de nos conectados
por arestas same-level e substituida por um unico no com:

```text
altitude preservada
proper parts agregados
filhos externos reanexados ao representante
```

Esse quociente preserva todos os conjuntos de limiar, porque nos de mesma
altitude aparecem exatamente nos mesmos limiares. Portanto, se `T~+` ja
reconstroi `I+`, a contracao same-level produz a forma canonica dessa mesma
arvore de componentes.

## Relacao com os nos `C_v^+`

O paper descreve a recomposicao usando nos locais por nivel, como `C_v^+`.
A implementacao atual nem sempre materializa explicitamente todos esses nos.
Em vez disso, a varredura por niveis reaproveita candidatos existentes, cria o
no de destino `b` quando necessario e, ao final, contrai arestas redundantes
de mesmo nivel.

Os surveys atuais indicam que:

```text
recomposicao local + criacao opcional de b + contracao same-level
```

tem produzido a mesma arvore canonica que o rebuild nos casos testados. Isso
sugere que parte dos `C_v^+` pode estar sendo representada implicitamente por
candidatos locais seguidos de contracao canonica.

O ponto teorico foi refinado em `docs/problem1_intermediate_gap_lemma.md`: para
o workspace single-`nodeCa`, se `v != b` e o bucket por valuation
`MergeNodes(v,-)` esta vazio, entao `C_v^+` e implicito. As arestas `Hasse+`
podem ser escritas diretamente sobre `Intervalle(C)`, usando o representante
canonico do proximo nivel materializado.

A ponte com a mutacao concreta esta em
`docs/problem1_explicit_region_mutation_proof.md`, que localiza essa
canonicalizacao no fechamento de `updateTreeWithExplicitRegion`.

A nota `docs/problem1_intermediate_level_gap_note.md` registra a
instrumentacao adicionada para medir essa diferenca. O resultado inicial e que
lacunas intermediarias sao frequentes mesmo quando a arvore final bate
estruturalmente o rebuild. A classificacao contra a arvore rebuild mostrou,
nos marcos testados, que essas lacunas intermediarias nao aparecem como nos
canonicos mesmo apos materializacao local de owners parciais; elas sao
implicitas e a canonicalizacao final remove a necessidade de materializa-las.

## Evidencia experimental atual

Com o contador `same_level_edges_collapsed`, observamos que a canonicalizacao
same-level e frequente:

```text
sintetico 8x8: 18995 contracoes, 6780/6780 eventos estruturais corretos
ICDAR 100x2 region16: 6936887 contracoes, 762954/762954 eventos estruturais corretos
ICDAR 64x64 controlado: 71799 contracoes, 1017/1017 eventos estruturais corretos
```

Assim, a proxima direcao teorica e formular a recomposicao local como:

```text
construir uma arvore possivelmente nao canonica localmente equivalente
aplicar contracao canonica de arestas de mesmo nivel
comparar com a arvore canonica de I+
```
