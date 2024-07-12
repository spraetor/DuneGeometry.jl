using QuadRules
using Documenter

makedocs(;
  modules=[QuadRules],
  author="Simon Praetorius",
  repo="https://github.com/spraetor/quadrules",
  format=Documenter.HTML(;
    canonical="https://spraetor.github.io/quadrules",
    assets=String[],
  ),
  pages=[
    "Home" => "index.md",
    "Library" => "Library.md",
  ],
)
