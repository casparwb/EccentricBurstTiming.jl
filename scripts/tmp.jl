function wolfram_to_julia(filename)

   expression = readlines(filename)
   expression = expression[2:end]
   expression[1] = replace(expression[1], "InputForm[" => "")
   expression[end] = expression[end][1:end-1]

   # filters = ["Pi" => "π",
   #            raw"""\[Omega]""" => "ω",
   #            raw"""\[CapitalOmega]""" => "Ω",
   #            raw"""\[Iota]""" => "ι",
   #            "Subscript[ω, -1 + i]" => "ωᵢ₋₁",
   #            "Subscript[Ω, -1 + i]" => "Ωᵢ₋₁",
   #            "Subscript[ι, -1 + i]" => "ιᵢ₋₁",
   #            "Subscript[V, 3]" => "V₃",
   #            "Subscript[m, 3]" => "m₃",
   #            "Subscript[e, -1 + i]" => "eᵢ₋₁",
   #            "Subscript[p, -1 + i]" => "pᵢ₋₁",
   #            "Subscript[ωᵢ₋₁, 3]" => "ω₃",
   #            "V₃" => "V₃ᵢ₋₁",
   #            "Sin" => "sin",
   #            "Cos" => "cos",
   #            "Sqrt" => "sqrt",
   #            "[" => "(", 
   #            "]" => ")",
   #            "M" => "m₁₂"
   #            ]

   filters = ["Pi" => "π",
              raw"""\[Omega]""" => "ω",
              raw"""\[CapitalOmega]""" => "Ω",
              raw"""\[Iota]""" => "ι",
              "Subscript[ω, 3]" => "ω₃",
              "Subscript[V, 3]" => "V₃",
              "Subscript[m, 3]" => "m₃",
              "Subscript[ωᵢ₋₁, 3]" => "ω₃",
              "Sin" => "sin",
              "Cos" => "cos",
              "Sqrt" => "sqrt",
              "[" => "(", 
              "]" => ")",
              "M" => "m₁₂"
              ]

   for i in eachindex(expression)
      e = expression[i]
      for fl in filters
         e = replace(e, fl)
      end
      expression[i] = e
   end

   expression = join(strip.(expression))
   write(replace(filename, ".wl" => "-jl.jl"), expression)
   nothing
end