// Resolves a DOI for every entry of scripts/data/tao-zhou-publications.json by
// querying Crossref and keeping only matches that agree on title, year and an author
// surname. Run with `node scripts/enrich-publication-dois.mjs` and commit the result;
// the site build never touches the network.
import fs from "node:fs"
import path from "node:path"
import { fileURLToPath } from "node:url"

const scriptDir = path.dirname(fileURLToPath(import.meta.url))
const dataPath = path.join(scriptDir, "data", "tao-zhou-publications.json")
const mailto = "freezeng-knowledge-base"

const normalize = (value) =>
  String(value)
    .toLowerCase()
    .replace(/[\u2010-\u2015]/g, "-")
    .replace(/[^a-z0-9]+/g, " ")
    .trim()

const tokens = (value) => new Set(normalize(value).split(" ").filter(Boolean))

function similarity(a, b) {
  const left = tokens(a)
  const right = tokens(b)
  if (left.size === 0 || right.size === 0) return 0
  let shared = 0
  for (const token of left) if (right.has(token)) shared += 1
  return (2 * shared) / (left.size + right.size)
}

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms))

async function crossref(paper) {
  const url = new URL("https://api.crossref.org/works")
  url.searchParams.set("query.bibliographic", paper.title)
  url.searchParams.set("rows", "5")
  url.searchParams.set("select", "title,DOI,container-title,issued,author")
  url.searchParams.set("mailto", mailto)
  const response = await fetch(url, { headers: { "User-Agent": `freezeng-kb (${mailto})` } })
  if (!response.ok) throw new Error(`Crossref returned HTTP ${response.status}`)
  const payload = await response.json()
  return payload.message?.items ?? []
}

function score(paper, item) {
  const title = Array.isArray(item.title) ? item.title[0] : item.title
  const titleScore = similarity(paper.title, title ?? "")
  const itemYear = item.issued?.["date-parts"]?.[0]?.[0]
  const yearPenalty = itemYear && Math.abs(itemYear - paper.year) > 2 ? 0.25 : 0
  const surnames = new Set(
    (item.author ?? []).map((author) => normalize(author.family ?? "")).filter(Boolean),
  )
  const wanted = paper.authors
    .split(",")
    .map((name) => normalize(name.trim().split(/\s+/).pop() ?? ""))
    .filter(Boolean)
  const authorHits = wanted.filter((surname) => surnames.has(surname)).length
  const authorScore = wanted.length ? authorHits / wanted.length : 0
  return { total: 0.72 * titleScore - yearPenalty + 0.28 * authorScore, titleScore, authorScore }
}

const data = JSON.parse(fs.readFileSync(dataPath, "utf8"))
const unresolved = []

for (const paper of data.papers) {
  if (paper.doi) continue
  if (paper.status === "submitted" && !paper.detail.startsWith("arXiv")) {
    unresolved.push(`${paper.n} (preprint, no DOI expected)`)
    continue
  }
  try {
    const items = await crossref(paper)
    let best = null
    for (const item of items) {
      const candidate = score(paper, item)
      if (!best || candidate.total > best.score.total) best = { item, score: candidate }
    }
    if (best && best.score.titleScore >= 0.7 && best.score.authorScore >= 0.5) {
      paper.doi = best.item.DOI
    } else {
      unresolved.push(`${paper.n} ${paper.title}`)
    }
  } catch (error) {
    unresolved.push(`${paper.n} ${paper.title} (${error.message})`)
  }
  await sleep(2500)
}

fs.writeFileSync(dataPath, `${JSON.stringify(data, null, 2)}\n`, "utf8")
const resolved = data.papers.filter((paper) => paper.doi).length
console.log(`Resolved ${resolved}/${data.papers.length} DOIs.`)
if (unresolved.length) console.log(`Unresolved:\n  ${unresolved.join("\n  ")}`)
