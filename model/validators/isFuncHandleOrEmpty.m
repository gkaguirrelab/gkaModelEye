function isFuncHandleOrEmpty(f)
    isFcn = isa(f, 'function_handle');
    if ~isFcn && ~isempty(f)
        error("Argument must be a function handle or empty")
    end
end