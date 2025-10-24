function wolfram_to_julia(e)

   filters = ["Pi" => π,
              raw"""\[Omega]""" => "ω",
              raw"""\[CapitalOmega]""" => "Ω",
              raw"""\[Iota]""" => "ι",
              "Subscript[ω, -1 + i]" => "ωᵢ₋₁",
              "Subscript[Ω, -1 + i]" => "Ωᵢ₋₁",
              "Subscript[ι, -1 + i]" => "ιᵢ₋₁",
              "Subscript[V, 3]" => "V₃",
              "Subscript[m, 3]" => "m₃",
              "Subscript[e, -1 + i]" => "eᵢ₋₁",
              "Subscript[p, -1 + i]" => "pᵢ₋₁",
              "Subscript[ωᵢ₋₁, 3]" => "ω₃",
              "Sin" => "sin",
              "Cos" => "cos",
              "[" => "(", 
              "]" => ")"
              ]

   for fl in filters
      e = replace(e, fl)
   end
   e
end