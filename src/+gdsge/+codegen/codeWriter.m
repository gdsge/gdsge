classdef codeWriter < handle
    % CODEWRITER  Minimal indented line buffer for MATLAB code emission.
    %   add(fmt, ...)  sprintf-formatted line at the current indent (4 sp/level)
    %   addRaw(text)   verbatim multi-line block, no indentation applied
    %   blank()        empty line
    %   in()/out()     indent level +/- 1
    %   str()          all lines joined with \n
    properties (Access = private)
        lines = {};
        level = 0;
    end
    methods
        function add(w, fmt, varargin)
            w.lines{end+1} = [repmat(' ', 1, 4*w.level) sprintf(fmt, varargin{:})];
        end
        function addRaw(w, text)
            text = strrep(text, sprintf('\r\n'), newline);
            parts = strsplit(text, newline, 'CollapseDelimiters', false);
            w.lines = [w.lines, parts];
        end
        function blank(w)
            w.lines{end+1} = '';
        end
        function in(w),  w.level = w.level + 1; end
        function out(w), w.level = max(0, w.level - 1); end
        function s = str(w)
            s = strjoin(w.lines, newline);
        end
    end
end
