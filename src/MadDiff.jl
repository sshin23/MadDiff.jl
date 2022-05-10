module MadDiff

using Reexport: @reexport

@reexport using MadDiffCore
@reexport using MadDiffModels
@reexport using MadDiffMOI

using MadDiffModels: Model
using MadDiffMOI: MadDiffAutomaticDifferentiation

export MadDiffAutomaticDifferentiation

end
