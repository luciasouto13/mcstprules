# Plot the Core Region and Cost Allocations

This function generates a geometric visualization of the *core* region
and the imputation set for cooperative games with 2, 3, or 4 players. It
also allows plotting multiple cost allocation rules simultaneously to
visually check their stability.

## Usage

``` r
mcstCorePlot(game = NULL, v = NULL, allocations = list(), titles = TRUE)
```

## Arguments

- game:

  list; an object containing a cooperative game (e.g., from
  [`mcstGamePrivate`](https://luciasouto13.github.io/mcstprules/reference/mcstGamePrivate.md)).
  If provided and `allocations` is empty, its allocation and solution
  concept are used by default. Default is `NULL`.

- v:

  numeric vector; the characteristic function values \\v(S)\\ for all
  \\2^n - 1\\ coalitions, ordered by cardinality. Required if `game` is
  `NULL`.

- allocations:

  list; a named list of numeric vectors representing different cost
  allocations to be plotted as points. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- titles:

  logical; if `TRUE` (default), adds a main title and subtitle/legend
  information to the plot.

## Value

For \\n = 2\\ or \\n = 3\\, the function draws a base R plot. For \\n =
4\\, it returns an interactive `plotly` object containing the 3D
visualization.

## Details

The geometric representation of the *core* and the imputation set
depends on the number of players \\n\\:

1.  \\n = 2\\ players: the imputation set is represented as a 1D line
    segment where \\x_1 + x_2 = v(N)\\. The *core* is highlighted as a
    sub-segment within it, delimited by individual rationality.

2.  \\n = 3\\ players: allocations are projected into 2D barycentric
    coordinates (simplex). The imputation set forms a large triangle,
    and the *core* is plotted inside as a bounded convex polygon.

3.  \\n = 4\\ players: the function switches to an interactive 3D plot
    using the `plotly` package. The imputation set is drawn as a
    transparent 3D tetrahedron (mesh), and the *core* is represented as
    a solid 3D polytope inside it.

If an allocation lies inside the green *core* region, it means the
allocation is coalitionally stable (no group of players has incentives
to defect).

The function automatically extracts the default allocation if a `game`
object is passed. Alternatively, a custom named `list` of allocations
can be provided to compare different rules (e.g., Bird, Shapley value,
folk) on the same graph.

## See also

[`mcstCoreCheck`](https://luciasouto13.github.io/mcstprules/reference/mcstCoreCheck.md)
for checking stability numbers and finding blocking coalitions.

[`mcstRules`](https://luciasouto13.github.io/mcstprules/reference/mcstRules.md)
for an overview of the available rules and analysis tools in the
package.

## Examples

``` r
# Example 1
vvals1 <- c(10, 12, 15, 20, 22, 25, 30) # 3 players: v(1), v(2), v(3), v(12), ...
allocs1 <- list(RuleA = c(8, 10, 12), RuleB = c(9, 10, 11), RuleC = c(7, 11, 12))
mcstCorePlot(v = vvals1, allocations = allocs1)

# \donttest{
# Example 2 (interactive 3D plot)
vvals2 <- c(10, 10, 10, 10, 18, 18, 18, 18, 18,
            18, 25, 25, 25, 25, 32) # 4 players: v(1), v(2), v(3), v(4), v(12), ...
allocs2 <- list(RuleA = c(8, 8, 8, 8), RuleB = c(10, 6, 11, 5))
mcstCorePlot(v = vvals2, allocations = allocs2)

{"x":{"visdat":{"1b521544c6f2":["function () ","plotlyVisDat"]},"cur_data":"1b521544c6f2","attrs":{"1b521544c6f2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"mesh3d","x":[10,10,10,2],"y":[10,10,2,10],"z":[10,2,10,10],"i":[2,1,1,1],"j":[0,0,2,2],"k":[3,3,3,0],"intensity":0,"showscale":false,"colorscale":[[0,"#E5E5E5"],[1,"#E5E5E5"]],"opacity":0.14999999999999999,"showlegend":false,"inherit":true},"1b521544c6f2.1":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"mesh3d","x":[8,10,10,8,10,8,7,7,7,7,7,7],"y":[7,7,7,7,8,10,8,10,10,8,7,7],"z":[7,7,8,10,7,7,7,7,8,10,10,8],"i":[7,0,4,9,10,10,10,10,4,4,4,4,9,9,9,9,9,9,9,9],"j":[8,11,2,10,0,2,2,0,0,7,7,0,2,8,4,4,7,11,7,10],"k":[5,6,1,3,1,1,3,11,1,6,5,6,3,5,5,2,6,6,8,11],"intensity":0,"showscale":false,"colorscale":[[0,"#74C476"],[1,"#74C476"]],"opacity":0.29999999999999999,"showlegend":false,"inherit":true},"1b521544c6f2.2":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":[8,10,10,8,10,8,7,7,7,7,7,7],"y":[7,7,7,7,8,10,8,10,10,8,7,7],"z":[7,7,8,10,7,7,7,7,8,10,10,8],"type":"scatter3d","mode":"markers","marker":{"size":1.2,"color":"black","opacity":0.59999999999999998},"showlegend":false,"inherit":true},"1b521544c6f2.3":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":8,"y":8,"z":8,"type":"scatter3d","mode":"markers","marker":{"size":2.5,"color":"#E41A1C","symbol":"square"},"name":"RuleA","inherit":true},"1b521544c6f2.4":{"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"x":10,"y":6,"z":11,"type":"scatter3d","mode":"markers","marker":{"size":2.5,"color":"#377EB8","symbol":"triangle-up"},"name":"RuleB","inherit":true}},"layout":{"margin":{"b":40,"l":60,"t":80,"r":100},"title":{"text":"<b>Core Region & Allocations<\/b><br><span style='font-size:12px'><span style='color:rgba(128,128,128,0.55);'>Imputation set<\/span> | <span style='color:rgba(49,163,84,0.55);'>Core<\/span><\/span>","x":0.5,"y":0.90000000000000002,"xanchor":"center"},"legend":{"y":0.5,"yanchor":"middle","x":1.05,"font":{"size":11},"itemwidth":15,"itemsizing":"constant","tracegroupgap":0},"scene":{"xaxis":{"title":"P1"},"yaxis":{"title":"P2"},"zaxis":{"title":"P3"},"aspectmode":"data"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"colorbar":{"title":"","ticklen":2},"colorscale":[[0,"#E5E5E5"],[1,"#E5E5E5"]],"showscale":false,"type":"mesh3d","x":[10,10,10,2],"y":[10,10,2,10],"z":[10,2,10,10],"i":[2,1,1,1],"j":[0,0,2,2],"k":[3,3,3,0],"intensity":[0],"opacity":0.14999999999999999,"showlegend":false,"frame":null},{"colorbar":{"title":"","ticklen":2},"colorscale":[[0,"#74C476"],[1,"#74C476"]],"showscale":false,"type":"mesh3d","x":[8,10,10,8,10,8,7,7,7,7,7,7],"y":[7,7,7,7,8,10,8,10,10,8,7,7],"z":[7,7,8,10,7,7,7,7,8,10,10,8],"i":[7,0,4,9,10,10,10,10,4,4,4,4,9,9,9,9,9,9,9,9],"j":[8,11,2,10,0,2,2,0,0,7,7,0,2,8,4,4,7,11,7,10],"k":[5,6,1,3,1,1,3,11,1,6,5,6,3,5,5,2,6,6,8,11],"intensity":[0],"opacity":0.29999999999999999,"showlegend":false,"frame":null},{"x":[8,10,10,8,10,8,7,7,7,7,7,7],"y":[7,7,7,7,8,10,8,10,10,8,7,7],"z":[7,7,8,10,7,7,7,7,8,10,10,8],"type":"scatter3d","mode":"markers","marker":{"color":"black","size":1.2,"opacity":0.59999999999999998,"line":{"color":"rgba(44,160,44,1)"}},"showlegend":false,"error_y":{"color":"rgba(44,160,44,1)"},"error_x":{"color":"rgba(44,160,44,1)"},"line":{"color":"rgba(44,160,44,1)"},"frame":null},{"x":[8],"y":[8],"z":[8],"type":"scatter3d","mode":"markers","marker":{"color":"#E41A1C","size":2.5,"symbol":"square","line":{"color":"rgba(214,39,40,1)"}},"name":"RuleA","error_y":{"color":"rgba(214,39,40,1)"},"error_x":{"color":"rgba(214,39,40,1)"},"line":{"color":"rgba(214,39,40,1)"},"frame":null},{"x":[10],"y":[6],"z":[11],"type":"scatter3d","mode":"markers","marker":{"color":"#377EB8","size":2.5,"symbol":"triangle-up","line":{"color":"rgba(148,103,189,1)"}},"name":"RuleB","error_y":{"color":"rgba(148,103,189,1)"},"error_x":{"color":"rgba(148,103,189,1)"},"line":{"color":"rgba(148,103,189,1)"},"frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}# }
```
