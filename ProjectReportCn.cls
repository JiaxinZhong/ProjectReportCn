\ProvidesClass{ProjectReportCn}[2024-8-29]


% -- 处理用户传入的语言选项 --
\newif\if@enlang % 新建一个布尔变量

\DeclareOption{en}{\@enlangtrue}
\DeclareOption{cn}{\@enlangfalse}
\ExecuteOptions{cn} % 默认中文
\ProcessOptions\relax

% 定义 \lang，供后面判断使用
\if@enlang
  \newcommand{\lang}{en}
\else
  \newcommand{\lang}{cn}
\fi

\LoadClass[openany]{report}

%% set the gemetric info for the document
\usepackage{geometry}
\geometry{
    paper = a4paper, % Change to letterpaper for US letter
    inner = 2.54cm, % Inner margin
    outer = 2.54cm, % Outer margin
    % bindingoffset = 0cm, % Binding offset
    top = 2.54cm, % Top margin
    bottom = 2.54cm, % Bottom margin
}


\usepackage{ifthen}

%% figures
\usepackage{graphicx}
% \usepackage{subcaption}

%% table
\usepackage{booktabs}


%% math
\usepackage{physics}
\usepackage{upgreek}
\usepackage{mathtools}

%% Set the footer
% see https://tex.stackexchange.com/questions/227/how-can-i-add-page-of-on-my-document
\usepackage{lastpage}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{} % clear existing header/footer entries
\fancyhead[L]{\leftmark}
\fancyhead[R]{\rightmark}
% \renewcommand{\headrule}{}
\ifthenelse{\equal{\lang}{cn}}{
    \fancyfoot[C]{第~{\thepage}~页 | 共~{\pageref*{LastPage}}~页}
}{
    \fancyfoot[C]{Page~{\thepage} of {\pageref*{LastPage}}}
}
\fancyfoot[L]{\titlestr}
\fancyfoot[R]{\datestr}
% 设置页脚线
% \renewcommand{\footrulewidth}{0.4pt}% default is 0pt

% note that 
%   \thispagestyle{firststyle}
% must be used below \maketitle
\fancypagestyle{firststyle}
{
    \fancyhf{}
    \ifthenelse{\equal{\lang}{cn}}{
        \fancyfoot[C]{第~{\thepage}~页 | 共~{\pageref*{LastPage}}~页}
    }{
        \fancyfoot[C]{Page~{\thepage} of {\pageref*{LastPage}}}
    }
    \fancyfoot[L]{\titlestr}
    \fancyfoot[R]{\datestr}
}
\fancypagestyle{plain}{%
    \fancyhf{}%
    \renewcommand{\headrule}{}
    \ifthenelse{\equal{\lang}{cn}}{
        \fancyfoot[C]{第~{\thepage}~页 | 共~{\pageref*{LastPage}}~页}
    }{
        \fancyfoot[C]{Page~{\thepage} of {\pageref*{LastPage}}}
    }
    \fancyfoot[L]{\titlestr}
    \fancyfoot[R]{\datestr}
}


% 设置 linspacing
\usepackage{setspace}
\onehalfspacing

% open book mark at arbitrary level
\usepackage[open, numbered]{bookmark}
% opt: numbered --- show numbers in the bookmark

% for some dummy text
\usepackage{lipsum}
\usepackage{zhlipsum}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable Chinese
\usepackage[UTF8]{ctex}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\usepackage{fontawesome5}
\usepackage{njuvisual}


%% font
\usepackage{metalogo}


% set heading styles
\usepackage{titlesec}
\titleformat{\chapter}{\centering\huge\bfseries}{第\,\thechapter\,章}{1em}{}

% 目录显示点
\usepackage[titles]{tocloft}
\renewcommand{\cftpartleader}{\cftdotfill{\cftdotsep}}
\renewcommand{\cftchapleader}{\cftdotfill{\cftdotsep}}
\renewcommand{\cftsecleader}{\cftdotfill{\cftdotsep}}
% 调整图索引字号
% \renewcommand{\cftfigfont}{\small}       % 条目文字（如 "图 1. 标题"）
% \renewcommand{\cftfigpagefont}{\small}  % 页码文字

% 修改目录
\usepackage{titletoc}

% see https://tex.stackexchange.com/questions/226481/appendix-section-title
\usepackage[title, titletoc]{appendix}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 插入代码块

% accsupp 宏包提供的 \BeginAccSupp/\EndAccSupp 来设定PDF内容的 ActualText
\usepackage{accsupp}
\newcommand{\emptyaccsupp}[1]{\BeginAccSupp{ActualText={}}#1\EndAccSupp{}}

\usepackage{listings}

\usepackage[numbered,framed]{matlab-prettifier}

% 默认字体改成 tt
\lstset{basicstyle=\ttfamily}

\newcommand{\codeinline}[1]{\texttt{\small#1}}


% 令 \lstinline!code!能分行
% https://stackoverflow.com/questions/1524289/verb-lstinline-text-flows-out-of-the-line-margins-without-line-break-in-latex
\lstset{
  literate={\_}{}{0\discretionary{\_}{}{\_}}%
}
% 自定义 MATLAB 的代码样式
\lstdefinestyle{matlab-custom}{
    % numbers=left,
    % stepnumber=1,
    % numbersep=10pt,
    tabsize=4,
    % showspaces=false,
    % showstringspaces=false,
    style=Matlab-editor,
    basicstyle=\ttfamily\small,
    belowcaptionskip=8pt, % 增大caption与代码之间的间距
    numberstyle=\tiny\emptyaccsupp, % 将行号包装为 ActualText 为空。关键作用：PDF中代码复制时不会包含行号
}

% Python 自定义代码风格
\lstdefinestyle{python-custom}{
    belowcaptionskip=8pt, % 增大caption与代码之间的间距
    language=Python,                                % 设置语言为 Python
    tabsize=4,                                      % 制表符宽度
    basicstyle=\ttfamily\small,                     % 设置等宽字体及字号，与 matlab-custom 保持一致
    keywordstyle=\color{blue!80!black}\bfseries,    % Python 保留关键字（如 def, for）高亮为蓝色粗体
    commentstyle=\color{gray!70!black}\itshape,     % 注释为灰色斜体
    stringstyle=\color{red!70!black},               % 字符串为略暗的红色
    identifierstyle=\color{black},                  % 变量名为黑色（默认）
    emph={np, plt, scipy, arange, linspace, array, exp, log, sin, cos, tan, abs},
    emphstyle=\color{orange!80!black},              % 强调常用科学函数，橙色显示
    showstringspaces=false,                         % 不特殊显示空格
    breaklines=true,                                % 自动换行，适配长表达式
    breakatwhitespace=true,                         % 仅在空白处分行，避免打断变量名
    numbers=left,                                   % 行号放在左侧
    numberstyle=\tiny\emptyaccsupp,                 % 行号为小字号，并避免复制时出现
    frame=single,                                   % 单线框包围代码
    framerule=0.5pt,                                % 边框线宽
    rulecolor=\color{gray!40},                      % 边框为浅灰色
    backgroundcolor=\color{gray!2},                 % 代码块背景为非常浅的灰色
    captionpos=t,                                   % caption（标题）显示在顶部
    xleftmargin=0em,                                % 左边距，用于对齐编号和代码体
    xrightmargin=0em                                % 右边距，防止代码太贴边
}


\ifthenelse{\equal{\lang}{cn}}{
    \renewcommand{\lstlistingname}{{\color{njupurple}\faFileCode[regular]} 代码}% Listing -> 代码
}{
    \renewcommand{\lstlistingname}{{\color{njupurple}\faFileCode[regular]} Code}% Listing -> 代码
}

% 将代码的字体设置为 Consolas
\setmonofont{Consolas}

% 自定义代码的标题
\newcommand{\codecaptionwithlink}[2]{
    #1\\
    \href{run:#2}{\color{njupurple}\faFile*[regular]}~\codeinline{#2}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 设置表格的行间距为1.5倍
% https://www.overleaf.com/learn/latex/Questions/How_do_I_change_column_or_row_separation_in_LaTeX_tables%3F
\renewcommand{\arraystretch}{1.5}


% listings的caption居中
\usepackage[justification=centering]{caption}

% 比较好看的字体
\usepackage[OT2,T1]{fontenc}
\DeclareSymbolFont{cyrletters}{OT2}{wncyr}{m}{n}
\DeclareMathSymbol{\Sha}{\mathalpha}{cyrletters}{"58}
% \usepackage{tgtermes}
\usepackage{fontspec} % 用于 XeLaTeX 或 LuaLaTeX
\setmainfont{TeX Gyre Termes} % 调用系统中可用的 TeX Gyre Termes 字体

\usepackage{xcolor}

% 支持 itemize 环境下的 itemsep
\usepackage{enumitem}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usepackage[acronym, shortcuts]{glossaries} % 缩写索引功能
\setglossarystyle{list} % 带页码的列表样式
\setacronymstyle{short-long}
\makeglossaries % 激活索引生成
\renewcommand*{\glstextformat}[1]{\textcolor{black}{#1}} % 把颜色改成黑色
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 打开书签
\usepackage[open, numbered]{bookmark}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 按键风格1
\usepackage[most]{tcolorbox}
\tcbuselibrary{listingsutf8}  % 支持 UTF-8
% 自定义按键样式
\newtcbox{\key}{on line,
  colback=gray!20, colframe=black,
  boxrule=0.4pt, arc=3pt, outer arc=1pt,
  boxsep=1pt, left=3pt, right=3pt, top=2pt, bottom=2pt,
  fontupper=\ttfamily\bfseries
}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%