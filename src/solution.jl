import Base.show, Base.print

mutable struct Route
   departure::Int64
   vertices::Array{Int}
end

mutable struct Solution
   cost::Float64
   routes::Array{Route}
end


# build Solution from the variables x
function getsolution(data::DataVRPRD, x, objval, optimizer)
   A, dim = arcs(data), n(data)+1
   adj_list = [[] for i in 1:dim] 
   for a in A
      val = get_value(optimizer, x[a])
      if val > 0.5
         push!(adj_list[a[1]+1], a[2])
      end
   end

   visited, routes = [false for i in 2:dim], []
   for i in adj_list[1]
      if !visited[i]
         r, prev = [], 0
         push!(r, i)
         visited[i] = true
         length(adj_list[i+1]) != 1 && error("Problem trying to recover the route from the x values. "*
                                             "Customer $i has $(length(adj_list[i+1])) outcoming arcs.")
         next, prev = adj_list[i+1][1], i
         maxit, it = dim, 0
         while next != 0 && it < maxit
            length(adj_list[next+1]) != 1 && error("Problem trying to recover the route from the x values. "* 
                                                   "Customer $next has $(length(adj_list[next+1])) outcoming arcs.")
            push!(r, next)
            visited[next] = true
            aux = next
            next, prev = adj_list[next+1][1], aux
            it += 1
         end
         (it == maxit) && error("Problem trying to recover the route from the x values. "*
                                "Route cannot be recovered because the return to depot is never reached")
         push!(routes, r)
      end
   end 
   !isempty(filter(y->y==false,visited)) && error("Problem trying to recover the route from the x values. "*
                              "At least one vertex was not visited or there are subtours in the solution x.")

   all_release_dates = []
   for i in 1:n(data)
      if(! (l(data, i) in all_release_dates) )
         push!(all_release_dates, l(data, i))
      end
   end

   sort!(all_release_dates)

   new_routes = []
   for r in routes
      curr_route = Route(0,r)
      for release_date in all_release_dates
         flag = 1
         for i in r
            if(l(data, i) == release_date)
               flag = 1
            elseif (l(data, i) < release_date && (release_date + t(data, (i, 0)) <= u(data,i)))
               flag = 1
            else
               flag = 0
            end

            if flag == 0
               break
            end
         end
         if flag == 1
            curr_route.departure = release_date
            break
         end
      end
      push!(new_routes,curr_route)
   end


   sol = Solution(objval, new_routes)
   sort!(sol.routes,by = r -> r.departure)
   return sol
end

function check_solution(data::DataVRPRD,solution)
   
   for (i,r) in enumerate(solution.routes)
      invalid = 0
      v = 0
      cur_time = r.departure
      println("Route $i: depart at $(cur_time)")

      for w in r.vertices
         cur_time += t(data,(v,w))
         println("On $w at $(cur_time)")
         if(cur_time > u(data,w))
            invalid = 1
            break
         end

         v = w
      end

      cur_time += t(data,(v,0))
      println("On 0 at $(cur_time)")
      if(cur_time > u(data,0))
         invalid = 1
         break
      end

      print("Route $i: ")
      if invalid == 1
         println("INVALID\n")
         break
      else
         println("OK\n")
      end
   end
end

function print_routes(solution)

   for (i,r) in enumerate(solution.routes)
      print("Route #$i - Departure at $(r.departure): ")
      for j in r.vertices
         print("$j ")
      end
      println()
   end
end

contains(p, s) = findnext(s, p, 1) != nothing
# read solution from file
function readsolution(solpath)
   str = read(solpath, String)
   breaks_in = [' '; ':'; '\n';'\t';'\r']
   aux = split(str, breaks_in; limit=0, keepempty=false) 
   sol = Solution(0, [])
   j = 3
   while j <= length(aux)
      r = []
      while j <= length(aux)
         push!(r, parse(Int, aux[j]))
         j += 1
         if contains(lowercase(aux[j]), "cost") || contains(lowercase(aux[j]), "route")
            break
         end
      end
      push!(sol.routes, r)
      if contains(lowercase(aux[j]), "cost")
         sol.cost = parse(Float64, aux[j+1])
         return sol
      end
      j += 2 # skip "Route" and "#j:" elements
   end
   error("The solution file was not read successfully. This format is not recognized.")
   return sol
end

# write solution in a file
function writesolution(solpath, solution)
   open(solpath, "w") do f
      for (i,r) in enumerate(solution.routes)
         write(f, "Route #$i - Departure at $(r.departure): ")
         for j in r.vertices
            write(f, "$j ") 
         end
         write(f, "\n")
      end
      write(f, "Cost $(solution.cost)\n")
   end
end

# write solution as TikZ figure (.tex) 
function drawsolution(tikzpath, data, solution)
   open(tikzpath, "w") do f
      write(f,"\\documentclass[crop,tikz]{standalone}\n")
      write(f, "\\usetikzlibrary{arrows}\n\\usetikzlibrary{shapes.geometric}\n\\usetikzlibrary {positioning}"*
               "\n\\pgfdeclarelayer{back}\n\\pgfsetlayers{back,main}\n\\makeatletter\n\\pgfkeys{%"*
               "\n/tikz/on layer/.code={\n\\def\\tikz@path@do@at@end{\\endpgfonlayer\\endgroup\\tikz@path@do@at@end}%"*
               "\n\\pgfonlayer{#1}\\begingroup%\n }%\n}\n\\makeatother\n\\begin{document}\n")

      # get limits to draw
      pos_x_vals = [i.pos_x for i in data.G′.V′]
      pos_y_vals = [i.pos_y for i in data.G′.V′]
      scale_fac = 1/(max(maximum(pos_x_vals),maximum(pos_y_vals))/10)
      write(f,"\\begin{tikzpicture}[thick, scale=1, every node/.style={minimum size=0.6cm, scale=0.4}, triangle/.style = {fill=white, regular polygon, regular polygon sides=6, scale=1.1, inner sep=0cm}]\n")
      for i in data.G′.V′
         x_plot, y_plot = scale_fac*i.pos_x, scale_fac*i.pos_y
         if i.id_vertex == 0 # plot depot
            write(f, "\t\\node[draw, line width=0.1mm, rectangle, fill=yellow, inner sep=0.05cm, scale=0.9] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
         elseif i.id_vertex <= n(data)
            write(f, "\t\\node[draw, line width=0.1mm, circle, fill=white, inner sep=0.05cm] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
         else
            write(f, "\t\\node[draw, line width=0.1mm, triangle, fill=white] (v$(i.id_vertex)) at ($(x_plot),$(y_plot)) {\\footnotesize $(i.id_vertex)};\n")
         end
      end
      write(f, "\\begin{pgfonlayer}{back}\n")


      for (idr,r) in enumerate(solution.routes)
         prev = 0
         for i in r.vertices
            a = (prev,i)
            edge_style = (prev == 0 || i == 0) ? "dashed,line width=0.2pt,opacity=.4" : "line width=0.4pt"
            write(f, "\t\\path[->,$(edge_style)] (v$(a[1])) edge node {} (v$(a[2]));\n")
            prev = i
         end
         write(f, "\t\\path[->,dashed,line width=0.2pt,opacity=.4] (v$prev) edge node {} (v0);\n")
      end
      write(f, "\\end{pgfonlayer}\n\\end{tikzpicture}\n")
      write(f, "\\end{document}\n")
   end
end
