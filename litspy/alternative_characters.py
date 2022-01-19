# EPMC can handle capitalisation, so it is not necessary to replace greek characters with all capitalisation variations
# however, many characters have e.g. greek and mathematical variations that look similar but are different
greek_dict = {"alpha": ["α", "𝛂", "𝛼"],
              "beta": ["β", "ϐ", "𝛽", "ᵝ"],
              "gamma": ["γ", "𝛄", "ℽ", "𝛾"],
              "delta": ["δ", "𝛿", "ẟ"],
              "epsilon": ["ε", "ɛ", "ϵ"],
              "zeta": ["ζ", "𝛇"],
              "eta": ["η"],
              "theta": ["Θ", "ϑ", "Ѳ"],  # Θ, ϴ, and θ are equivalent in EPMC
              "iota": ["Ι", "Ɩ"],
              "kappa": ["Κ", "ϰ"],
              "lambda": ["Λ"],
              "mu": ["Μ", "µ", "𝜇", "𝝁"],  # lowercase mu and the symbol for the micro sign look the same in many fonts
              "nu": ["Ν", "𝜈"],
              "xi": ["ξ"],
              "omicron": ["Ο"],
              "pi": ["Π", "ϖ", "𝜋"],
              "rho": ["Ρ"],
              "sigma": ["Σ", "ς", "𝜎"],
              "tau": ["Τ"],
              "upsilon": ["Υ", "ϒ"],
              "phi": ["φ", "ϕ", "Ф"],
              "chi": ["χ"],
              "psi": ["ψ", "𝛹"],
              "omega": ["Ω", "ѡ"]}

numerals = ["I", "X", "V"]

hyphens = ["-", "–", "—", "‑"]
