# Enable Powerlevel10k instant prompt. Should stay close to the top of ~/.zshrc.
# Initialization code that may require console input (password prompts, [y/n]
# confirmations, etc.) must go above this block; everything else may go below.
if [[ -r "${XDG_CACHE_HOME:-$HOME/.cache}/p10k-instant-prompt-${(%):-%n}.zsh" ]]; then
  source "${XDG_CACHE_HOME:-$HOME/.cache}/p10k-instant-prompt-${(%):-%n}.zsh"
fi

  #color{{{
  autoload colors
  colors
   
  for color in RED GREEN YELLOW BLUE MAGENTA CYAN WHITE; do
  eval _$color='%{$terminfo[bold]$fg[${(L)color}]%}'
  eval $color='%{$fg[${(L)color}]%}'
  (( count = $count + 1 ))
  done
  FINISH="%{$terminfo[sgr0]%}"
  #}}}

  # If you come from bash you might have to change your $PATH.
  # export PATH=$HOME/bin:/usr/local/bin:$PATH

  # Path to your oh-my-zsh installation.
  export ZSH="$HOME/.oh-my-zsh"

  # Set name of the theme to load --- if set to "random", it will
  # load a random theme each time oh-my-zsh is loaded, in which case,
  # to know which specific one was loaded, run: echo $RANDOM_THEME
  # See https://github.com/robbyrussell/oh-my-zsh/wiki/Themes
  ZSH_THEME="powerlevel10k/powerlevel10k"

  # Set list of themes to pick from when loading at random
  # Setting this variable when ZSH_THEME=random will cause zsh to load
  # a theme from this variable instead of looking in ~/.oh-my-zsh/themes/
  # If set to an empty array, this variable will have no effect.
  # ZSH_THEME_RANDOM_CANDIDATES=( "robbyrussell" "agnoster" )

  # Uncomment the following line to use case-sensitive completion.
  # CASE_SENSITIVE="true"

  # Uncomment the following line to use hyphen-insensitive completion.
  # Case-sensitive completion must be off. _ and - will be interchangeable.
  # HYPHEN_INSENSITIVE="true"

  # Uncomment the following line to disable bi-weekly auto-update checks.
  # DISABLE_AUTO_UPDATE="true"

  # Uncomment the following line to automatically update without prompting.
  # DISABLE_UPDATE_PROMPT="true"

  # Uncomment the following line to change how often to auto-update (in days).
  # export UPDATE_ZSH_DAYS=13

  # Uncomment the following line if pasting URLs and other text is messed up.
  # DISABLE_MAGIC_FUNCTIONS=true

  # Uncomment the following line to disable colors in ls.
  # DISABLE_LS_COLORS="true"

  # Uncomment the following line to disable auto-setting terminal title.
  # DISABLE_AUTO_TITLE="true"

  # Uncomment the following line to enable command auto-correction. 用The FUCK插件代替
  # ENABLE_CORRECTION="true"

  # Uncomment the following line to display red dots whilst waiting for completion.
  # COMPLETION_WAITING_DOTS="true"

  # Uncomment the following line if you want to disable marking untracked files
  # under VCS as dirty. This makes repository status check for large repositories
  # much, much faster.
  # DISABLE_UNTRACKED_FILES_DIRTY="true"

  # Uncomment the following line if you want to change the command execution time
  # stamp shown in the history command output.
  # You can set one of the optional three formats:
  # "mm/dd/yyyy"|"dd.mm.yyyy"|"yyyy-mm-dd"
  # or set a custom format using the strftime function format specifications,
  # see 'man strftime' for details.
  HIST_STAMPS="mm/dd/yyyy"

  # Would you like to use another custom folder than $ZSH/custom?
  # ZSH_CUSTOM=/path/to/new-custom-folder

  # Which plugins would you like to load?
  # Standard plugins can be found in ~/.oh-my-zsh/plugins/*
  # Custom plugins may be added to ~/.oh-my-zsh/custom/plugins/
  # Example format: plugins=(rails git textmate ruby lighthouse)
  # Add wisely, as too many plugins slow down shell startup.
  # z命令快速跳转目录     x命令解压一切文件
  plugins=(
    git
    z
    zsh-autosuggestions
    extract
    zsh-syntax-highlighting
    last-working-dir
  )

  source $ZSH/oh-my-zsh.sh

  # User configuration

  # export MANPATH="/usr/local/man:$MANPATH"

  # You may need to manually set your language environment
  # export LANG=en_US.UTF-8

  # Preferred editor for local and remote sessions
  # if [[ -n $SSH_CONNECTION ]]; then
  #   export EDITOR='vim'
  # else
  #   export EDITOR='mvim'
  # fi

  # Compilation flags
  # export ARCHFLAGS="-arch x86_64"

  # Set personal aliases, overriding those provided by oh-my-zsh libs,
  # plugins, and themes. Aliases can be placed here, though oh-my-zsh
  # users are encouraged to define aliases within the ZSH_CUSTOM folder.
  # For a full list of active aliases, run `alias`.
  #
  # Example aliases 别名
  # alias zshconfig="mate ~/.zshrc"
  # alias ohmyzsh="mate ~/.oh-my-zsh"
  alias check='nvidia-smi'
  alias ta='tmux attach'
  alias cp='cp -i'
  alias mv='mv -i'
  alias rm='rm -i'
  alias ls='ls -F --color=auto'
  alias ll='ls -al'
  alias grep='grep --color=auto'
  alias la='ls -a'
  # list of alias
  ## Colorize the ls output ##
  alias ls='ls --color=auto'

  ## Use a long listing format ##
  alias ll='ls -la'

  ## get rid of command not found ##
  alias cd..='cd ..'

  ## a quick way to get out of current directory ##
  alias ..='cd ..'
  alias ...='cd ../../../'
  alias ....='cd ../../../../'
  alias .....='cd ../../../../'
  alias .4='cd ../../../../'
  alias .5='cd ../../../../..'

  ## Colorize the grep command output for ease of use (good for log files)##
  alias grep='grep --color=auto'
  alias egrep='egrep --color=auto'
  alias fgrep='fgrep --color=auto'

  ## this one saved by butt so many times ##
  alias wget='wget -c'

  #路径别名 {{{
  #进入相应的路径时只要 cd ~xxx
  hash -d H="/home"
  #}}}

#   eval $(thefuck --alias)

  #命令提示符样式{{{
  RPROMPT=$(echo "$RED%D %T$FINISH")
  PROMPT=$(echo "$WHITE@SSH:  $CYAN%n@$YELLOW%M:$GREEN%/$_YELLOW>$FINISH ")
  #}}}

  #标题栏、任务栏样式{{{
  case $TERM in (*xterm*|*rxvt*|(dt|k|E)term)
  precmd () { print -Pn "\e]0;%n@%M//%/\a" }
  preexec () { print -Pn "\e]0;%n@%M//%/\ $1\a" }
  ;;
  esac
  #}}}
  
  #编辑器
  export EDITOR=vim

  #关于历史纪录的配置 {{{
  #历史纪录条目数量
  export HISTSIZE=10000
  #注销后保存的历史纪录条目数量
  export SAVEHIST=10000
  #历史纪录文件
  export HISTFILE=~/.zsh_history
  #以附加的方式写入历史纪录
  setopt INC_APPEND_HISTORY
  #如果连续输入的命令相同，历史纪录中只保留一个
  setopt HIST_IGNORE_DUPS
  #为历史纪录中的命令添加时间戳
  setopt EXTENDED_HISTORY      
   
  #启用 cd 命令的历史纪录，cd -[TAB]进入历史路径
  setopt AUTO_PUSHD
  #相同的历史路径只保留一个
  setopt PUSHD_IGNORE_DUPS
   
  #在命令前添加空格，不将此命令添加到纪录文件中
  #setopt HIST_IGNORE_SPACE
  #}}}

  ##空行(光标在行首)补全 "cd " {{{
  user-complete(){
  case $BUFFER in
  "" )                       # 空行填入 "cd "
  BUFFER="cd "
  zle end-of-line
  zle expand-or-complete
  ;;
  "cd --" )                  # "cd --" 替换为 "cd +"
  BUFFER="cd +"
  zle end-of-line
  zle expand-or-complete
  ;;
  "cd +-" )                  # "cd +-" 替换为 "cd -"
  BUFFER="cd -"
  zle end-of-line
  zle expand-or-complete
  ;;
  * )
  zle expand-or-complete
  ;;
  esac
  }
  zle -N user-complete
  bindkey "\t" user-complete
  #}}}

# To customize prompt, run `p10k configure` or edit ~/.p10k.zsh.
[[ ! -f ~/.p10k.zsh ]] || source ~/.p10k.zsh

alias	le='less -mNS'
alias	conbase='conda activate'
alias	conde='conda deactivate'
alias   conact='conda activate'

alias   Rno='R --no-restore'  #for R start without restore current .RData