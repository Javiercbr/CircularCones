function config = load_config(filename)

    if nargin < 1
        filename = 'config.json';
    end

    fid = fopen(filename, 'r');

    if fid == -1
        error('Could not open configuration file: %s', filename);
    end

    cleanup = onCleanup(@() fclose(fid));

    raw = fread(fid, inf, '*char')';
    config = jsondecode(raw);

end