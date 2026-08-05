// Regenerates the bilingual publication catalogue pages from
// scripts/data/tao-zhou-publications.json. Run `node scripts/generate-publication-catalog.mjs`
// after editing the dataset. Everything outside the generated markers is hand-written and
// preserved.
import fs from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const repoRoot = path.resolve(scriptDir, "..")
const data = JSON.parse(
  fs.readFileSync(path.join(scriptDir, "data", "tao-zhou-publications.json"), "utf8"),
)

const BEGIN = "<!-- generated:catalog:begin -->"
const END = "<!-- generated:catalog:end -->"

const targets = {
  zh: {
    file: path.join(repoRoot, "content", "computational-mathematics", "paper-notes", "catalog.md"),
    root: "computational-mathematics/paper-notes",
    headers: ["编号", "论文", "发表信息", "精读页"],
    submitted: "投稿或预印本",
    published: "已发表",
    countLabel: (n) => `${n} 篇`,
    noNote: "—",
    themeHeading: (theme) => `${theme.zh}`,
  },
  en: {
    file: path.join(
      repoRoot,
      "content",
      "en",
      "computational-mathematics",
      "paper-notes",
      "catalog.md",
    ),
    root: "en/computational-mathematics/paper-notes",
    headers: ["No.", "Paper", "Publication", "Close reading"],
    submitted: "Submitted or preprint",
    published: "Published",
    countLabel: (n) => `${n} papers`,
    noNote: "—",
    themeHeading: (theme) => `${theme.en}`,
  },
}

const escapeCell = (value) => String(value).replaceAll("|", "\\|")

function paperTitleCell(paper) {
  const title = escapeCell(paper.title)
  return paper.doi ? `[${title}](https://doi.org/${paper.doi})` : title
}

function venueCell(paper) {
  const parts = [paper.venue]
  if (paper.detail) parts.push(paper.detail)
  parts.push(String(paper.year))
  return escapeCell(parts.join(", "))
}

function noteCell(paper, target, pageTitles) {
  const slug = data.themes.find((theme) => theme.id === paper.theme)?.slug
  if (!slug || !paper.page) return target.noNote
  const key = `${slug}/${paper.page}`
  const label = pageTitles[key]
  if (!label) return target.noNote
  // The alias separator must be escaped: an unescaped pipe would be read as a
  // table column separator, splitting the cell and leaving the wikilink
  // unparsed. Quartz's wikilink regex accepts `\|`.
  return `[[${target.root}/${key}\\|${label}]]`
}

function buildTable(rows, headers) {
  const lines = [`| ${headers.join(" | ")} |`, `| ${headers.map(() => "---").join(" | ")} |`]
  for (const row of rows) lines.push(`| ${row.join(" | ")} |`)
  return lines.join("\n")
}

function collectPageTitles(target, lang) {
  const titles = {}
  const base = path.dirname(target.file)
  for (const theme of data.themes) {
    const dir = path.join(base, theme.slug)
    if (!fs.existsSync(dir)) continue
    for (const entry of fs.readdirSync(dir)) {
      if (!entry.endsWith(".md") || entry === "index.md") continue
      const raw = fs.readFileSync(path.join(dir, entry), "utf8")
      const match = raw.match(/^---\r?\n[\s\S]*?^title:\s*(.+?)\r?$/m)
      if (!match) continue
      titles[`${theme.slug}/${entry.replace(/\.md$/, "")}`] = match[1]
        .trim()
        .replace(/^["']|["']$/g, "")
    }
  }
  return titles
}

let written = 0
for (const [lang, target] of Object.entries(targets)) {
  if (!fs.existsSync(target.file)) {
    console.warn(`Skipping ${path.relative(repoRoot, target.file)}: file not found.`)
    continue
  }
  const pageTitles = collectPageTitles(target, lang)
  const sections = []
  for (const theme of data.themes) {
    const papers = data.papers.filter((paper) => paper.theme === theme.id).sort((a, b) => b.n - a.n)
    const published = papers.filter((paper) => paper.status === "published")
    const submitted = papers.filter((paper) => paper.status === "submitted")
    const heading = `### ${target.themeHeading(theme)} · ${target.countLabel(papers.length)}`
    const blocks = [heading]
    if (submitted.length) {
      blocks.push(`**${target.submitted}**`)
      blocks.push(
        buildTable(
          submitted.map((paper) => [
            paper.n,
            paperTitleCell(paper),
            venueCell(paper),
            noteCell(paper, target, pageTitles),
          ]),
          target.headers,
        ),
      )
    }
    blocks.push(`**${target.published}**`)
    blocks.push(
      buildTable(
        published.map((paper) => [
          paper.n,
          paperTitleCell(paper),
          venueCell(paper),
          noteCell(paper, target, pageTitles),
        ]),
        target.headers,
      ),
    )
    sections.push(blocks.join("\n\n"))
  }

  const generated = [BEGIN, "", sections.join("\n\n"), "", END].join("\n")
  const raw = fs.readFileSync(target.file, "utf8")
  const start = raw.indexOf(BEGIN)
  const end = raw.indexOf(END)
  if (start === -1 || end === -1) {
    console.warn(`Skipping ${path.relative(repoRoot, target.file)}: markers not found.`)
    continue
  }
  const next = `${raw.slice(0, start)}${generated}${raw.slice(end + END.length)}`
  fs.writeFileSync(target.file, next, "utf8")
  written += 1
}

console.log(`Rewrote the generated catalogue block in ${written} file(s).`)
