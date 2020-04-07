# Screen and Tmux Best Pratice

The program `screen` allows you to detach your session from your current login session.
So in case you get disconnected your screen session will stay alive.

!!! hint
    You have to reconnect to screen on the machine that you started it.
    We thus recommend starting it only on the login nodes and **not** on a compute node.

## Start and terminat a screen session

You start a new `screen` session by

```terminal
$ screen
```
When you are in a screen session you can terminate it with

```terminal
$ exit
```
so its gone then.

## Detach a screen session

If you want to detach your screen session press `Ctrl+a d`

## List screen sessions

To list all your screen sessions run

```terminal
$ screen -ls

There is a screen on:
	2441.pts-1.med0236	(Detached)
1 Socket in /var/run/screen/S-kbentel.
```

## Reattach screen session

To reattach a screen session run

```terminal
$ screen -r screen_session_id
```

If you do not know the `screen_session_id` you can get it with `screen -ls`, e.g. `2441.pts-1.med0236` in the example above. You do not have to type the whole `screen_session_id` only as much as is necessary to identify it uniquely. In case there is only one screen session detached it is enough to run `screen -r`

## Kill a detached screen session

Sometimes it is necessary to kill a detached screen session. This is done with the command

```terminal
$ screen -X -S screen_session_id quit
```

## Multiple windows in a screen session

It is possible to have multiple windows in a screen session. So suppose you are logged into a screen session, these are the relevant shortcuts

```
new win:           Ctrl+a c
next/previous win: Ctrl+a n/p
```

To terminate a window just enter

```terminal
$ exit
```

## Configuration file

Here is a sensible screen configuration.
Save it as `~/.screenrc`.

- **TODO** [screenrc](files/screenrc)

## Fix a broken screen session

In case your screen session doesn't write to the terminal correctly,
i.e. the formatting of the output is broken, you can fix it by typing
to the terminal:

```terminal
$ tput smam
```
