function results = initialize_results(options)
    results.min_set = options.starting_set;
    results.nambiguous_best = Inf;
    results.nmissed_best = Inf;
end