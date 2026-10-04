//! Shared confirmation selector: arrow-key (or j/k) choice between options,
//! Enter confirms, Esc/Ctrl-C aborts. Defaults to the safe option (first).
//! Non-TTY callers get an automatic `false` unless the caller passed a
//! script bypass flag.

use crossterm::event::{Event, KeyCode, KeyEventKind};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode};

/// Interactive choice between `options` (first = safe default). Returns the
/// chosen index, or `None` on abort. Never raw-modes a non-TTY.
pub fn choose(prompt: &str, options: &[&str]) -> anyhow::Result<Option<usize>> {
    use std::io::IsTerminal as _;
    if !std::io::stdin().is_terminal() {
        return Ok(None);
    }
    enable_raw_mode()?;
    let mut sel = 0usize;
    let mut out = std::io::stdout();
    use std::io::Write as _;
    let result = (|| -> anyhow::Result<Option<usize>> {
        loop {
            let mut line = String::from(prompt);
            line.push('\n');
            for (i, o) in options.iter().enumerate() {
                if i == sel {
                    line.push_str(&format!("  > {o}\n"));
                } else {
                    line.push_str(&format!("    {o}\n"));
                }
            }
            write!(out, "{line}")?;
            out.flush()?;
            let ev = crossterm::event::read()?;
            match ev {
                Event::Key(k) if k.kind == KeyEventKind::Press => match k.code {
                    KeyCode::Up | KeyCode::Char('k') => sel = sel.saturating_sub(1),
                    KeyCode::Down | KeyCode::Char('j') => sel = (sel + 1).min(options.len() - 1),
                    KeyCode::Enter => {
                        writeln!(out)?;
                        out.flush()?;
                        return Ok(Some(sel));
                    }
                    KeyCode::Esc => {
                        writeln!(out)?;
                        out.flush()?;
                        return Ok(None);
                    }
                    _ => {}
                },
                _ => {}
            }
            // redraw: move cursor up over the menu
            write!(out, "\x1b[{}A\x1b[0J", options.len() + 1)?;
            out.flush()?;
        }
    })();
    disable_raw_mode()?;
    result
}

/// True when the user picked the confirm option. `labels` = [safe, danger].
pub fn confirmed(prompt: &str, labels: (&str, &str), yes: bool) -> anyhow::Result<bool> {
    if yes {
        return Ok(true);
    }
    match choose(prompt, &[labels.0, labels.1])? {
        Some(1) => Ok(true),
        _ => Ok(false),
    }
}
