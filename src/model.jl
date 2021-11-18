function build_model(data::DataVRPRD, app)

   A = arcs(data) 
   V = [i for i in 1:n(data)] 
   Q = veh_capacity(data)

   # Formulation
   vrprd = VrpModel()
   @variable(vrprd.formulation, x[a in A], Int)
   @objective(vrprd.formulation, Min, sum(c(data,a) * x[a] for a in A))
   @constraint(vrprd.formulation, indeg[i in V], sum(x[a] for a in A if a[2] == i) == 1.0)

   #println(vrprd.formulation)

   function buildgraph(release_date::Int)
      v_source = v_sink = 0
     
      V1 = [0]
      for i in 1:n(data)
         if(l(data, i) == release_date)
            push!(V1,i)
         elseif (l(data, i) < release_date && (release_date + t(data, (i, 0)) <= u(data,i)))
            push!(V1,i)
         end   
      end

      #println(V1)
      
      L, U = 0, upperBoundNbVehicles(data) # multiplicity

      G = VrpGraph(vrprd, V1, v_source, v_sink, (L, U))

      if app["enable_cap_res"]
         cap_res_id = add_resource!(G, main = true)
      end
      time_res_id = add_resource!(G, main = true)
      

      for v in V1
         if app["enable_cap_res"]
            set_resource_bounds!(G, v, cap_res_id, 0, Q)
         end

         
         set_resource_bounds!(G, v, time_res_id, 0, u(data, v) - release_date)
      end

      for (i,j) in A
         if(i in V1 && j in V1)
            arc_id = add_arc!(G, i, j)
            add_arc_var_mapping!(G, arc_id, x[(i,j)])

            if app["enable_cap_res"]
               set_arc_consumption!(G, arc_id, cap_res_id, d(data, j))
            end
            set_arc_consumption!(G, arc_id, time_res_id, t(data, (i, j)))
         end
      end

      return G,V1
   end
   
   all_release_dates = []
   for i in 1:n(data)
      if(! (l(data, i) in all_release_dates) )
         push!(all_release_dates, l(data, i))
      end
   end

   #println(all_release_dates)

   graphs, packing_set_vertex = [], [[] for v in 1:n(data)]
   k = 1
   for rd in all_release_dates
      G,V1 = buildgraph(Int(rd))
      add_graph!(vrprd, G)   
      println(G)
      push!(graphs,G)
      #println(V1, packing_set_vertex)
      for v in V1
         if(v == 0)
            continue
         end
         push!(packing_set_vertex[v],k)
      end
      k+=1
   end

   #define_elementarity_sets_distance_matrix!(vrprd, G, [[c(data, (i, j)) for j in V] for i in V])

   set_vertex_packing_sets!(vrprd, [[(graphs[k],v) for k in packing_set_vertex[v]]  for v in 1:n(data)])
   add_capacity_cut_separator!(vrprd, [ ( [(graphs[k],v) for k in packing_set_vertex[v]], Float64(d(data, v))) for v in 1:n(data)], Float64(Q))

   set_branching_priority!(vrprd, "x", 1)

   return (vrprd, x)
end
