using Coverage
using Printf

coverage = process_folder("src")
covered_lines, total_lines = get_summary(coverage)
coverage_pct = total_lines == 0 ? 0.0 : 100 * covered_lines / total_lines

LCOV.writefile("lcov.info", coverage)

open("coverage-summary.md", "w") do io
    println(io, "## Coverage report")
    @printf(io, "- Overall line coverage: **%.2f%%** (%d/%d lines)\n", coverage_pct, covered_lines, total_lines)
    println(io, "- Scope: `src/`")
    println(io, "- Format: `lcov.info` generated via Coverage.jl")
    println(io)
    println(io, "_This comment is CI-generated and updates on each run._")
end

open(ENV["GITHUB_STEP_SUMMARY"], "a") do io
    println(io, "## Coverage")
    @printf(io, "- Overall line coverage: **%.2f%%** (%d/%d lines)\n", coverage_pct, covered_lines, total_lines)
    println(io, "- Scope: `src/`")
end
