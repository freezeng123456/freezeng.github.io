import { QuartzComponent, QuartzComponentConstructor, QuartzComponentProps } from "./types"
import { FullSlug, joinSegments, pathToRoot } from "../util/path"

type LanguageFrontmatter = {
  lang?: string
  translation?: string
}

const LanguageSwitcher: QuartzComponent = ({ fileData }: QuartzComponentProps) => {
  const frontmatter = (fileData.frontmatter ?? {}) as LanguageFrontmatter
  const slug = fileData.slug ?? "index"
  const isEnglish = frontmatter.lang === "en" || slug === "en" || slug.startsWith("en/")
  const target = frontmatter.translation ?? (isEnglish ? "index" : "en")
  const baseDir = pathToRoot(slug as FullSlug)
  const href = joinSegments(baseDir, target as FullSlug)
  const label = isEnglish ? "切换到中文版" : "Switch to English"

  return (
    <a
      class="language-switcher internal"
      href={href}
      aria-label={label}
      title={label}
      data-current-language={isEnglish ? "en" : "zh"}
    >
      <span class={isEnglish ? "language-option" : "language-option active"} lang="zh">
        中
      </span>
      <span class="language-divider" aria-hidden="true">
        /
      </span>
      <span class={isEnglish ? "language-option active" : "language-option"} lang="en">
        EN
      </span>
    </a>
  )
}

LanguageSwitcher.css = `
.language-switcher {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 0.2rem;
  min-width: 3.85rem;
  height: 2rem;
  padding: 0 0.45rem;
  border: 1px solid transparent;
  border-radius: 0.45rem;
  color: var(--darkgray);
  font-family: var(--bodyFont);
  font-size: 0.78rem;
  font-weight: 600;
  line-height: 1;
  text-decoration: none;
  white-space: nowrap;
  transition:
    color 0.2s ease,
    background-color 0.2s ease,
    border-color 0.2s ease;
}

.language-switcher:hover {
  color: var(--dark);
  background: var(--highlight);
  border-color: var(--lightgray);
}

.language-switcher:focus-visible {
  outline: 2px solid var(--secondary);
  outline-offset: 2px;
}

.language-switcher .language-option {
  opacity: 0.52;
}

.language-switcher .language-option.active {
  color: var(--secondary);
  opacity: 1;
}

.language-switcher .language-divider {
  color: var(--gray);
  font-weight: 400;
}

@media (max-width: 800px) {
  .language-switcher {
    min-width: 3.5rem;
    padding-inline: 0.35rem;
  }
}
`

LanguageSwitcher.afterDOMLoaded = `
const localizeKnowledgeBaseUi = () => {
  const slug = document.body.dataset.slug ?? ""
  const isEnglish = slug === "en" || slug.startsWith("en/")
  const basePath = document.body.dataset.basepath ?? ""

  document.documentElement.lang = isEnglish ? "en" : "zh"

  const siteTitle = document.querySelector(".page-title a")
  if (siteTitle) {
    siteTitle.textContent = isEnglish ? "Freezeng Knowledge Base" : "Freezeng 知识库"
    siteTitle.setAttribute("href", basePath + (isEnglish ? "/en/" : "/"))
  }

  const switcher = document.querySelector(".language-switcher")
  if (switcher) {
    const label = isEnglish ? "切换到中文版" : "Switch to English"
    switcher.setAttribute("aria-label", label)
    switcher.setAttribute("title", label)
  }

  const searchButton = document.querySelector(".search-button")
  const searchButtonText = searchButton?.querySelector("p")
  const searchInput = document.querySelector(".search-bar")
  if (searchButton) searchButton.setAttribute("aria-label", isEnglish ? "Search" : "搜索")
  if (searchButtonText) searchButtonText.textContent = isEnglish ? "Search" : "搜索"
  if (searchInput) {
    searchInput.setAttribute("aria-label", isEnglish ? "Search the knowledge base" : "搜索知识库")
    searchInput.setAttribute("placeholder", isEnglish ? "Search" : "搜索")
  }

  const labels = [
    [".explorer-toggle h2", isEnglish ? "Explore" : "探索"],
    [".graph > h3", isEnglish ? "Graph View" : "关系图谱"],
    [".toc-header h3", isEnglish ? "Table of Contents" : "目录"],
    [".backlinks > h3", isEnglish ? "Backlinks" : "反向链接"],
    [".recent-notes > h3", isEnglish ? "Recent Notes" : "最近更新"],
  ]
  for (const [selector, text] of labels) {
    const element = document.querySelector(selector)
    if (element) element.textContent = text
  }

  const darkmode = document.querySelector(".darkmode")
  const readermode = document.querySelector(".readermode")
  if (darkmode) {
    darkmode.setAttribute("aria-label", isEnglish ? "Toggle color theme" : "切换颜色主题")
  }
  if (readermode) {
    readermode.setAttribute("aria-label", isEnglish ? "Reader mode" : "阅读模式")
  }

  const filterLanguageLinks = (root) => {
    if (!root) return

    if (root.classList.contains("recent-notes")) {
      let visibleCount = 0
      root.querySelectorAll(".recent-li").forEach((item) => {
        const link = item.querySelector("h3 a")
        if (!link) return
        const path = new URL(link.href, window.location.href).pathname
        const englishTarget =
          path === "/en" ||
          path.startsWith("/en/") ||
          path === basePath + "/en" ||
          path.startsWith(basePath + "/en/")
        const sameLanguage = isEnglish ? englishTarget : !englishTarget
        const visible = sameLanguage && visibleCount < 5
        item.toggleAttribute("hidden", !visible)
        if (visible) visibleCount += 1
      })
      return
    }

    const resultCards = Array.from(root.querySelectorAll("a.result-card"))
    const matchesLanguage = (link) => {
      const path = new URL(link.href, window.location.href).pathname
      const englishTarget =
        path === "/en" ||
        path.startsWith("/en/") ||
        path === basePath + "/en" ||
        path.startsWith(basePath + "/en/")
      return isEnglish ? englishTarget : !englishTarget
    }
    const hasSameLanguageResult = resultCards.some(matchesLanguage)
    resultCards.forEach((link) => {
      link.toggleAttribute("hidden", hasSameLanguageResult && !matchesLanguage(link))
    })
  }

  const searchLayout = document.querySelector(".search-layout")
  const recentNotes = document.querySelector(".recent-notes")
  filterLanguageLinks(searchLayout)
  filterLanguageLinks(recentNotes)

  const orderTimeParallelChapters = (list, linkSelector) => {
    if (!list) return
    const chapterItems = Array.from(list.children)
      .map((item) => {
        if (!(item instanceof HTMLElement)) return null
        const link = item.querySelector(linkSelector)
        if (!(link instanceof HTMLAnchorElement)) return null
        const match = new URL(link.href, window.location.href).pathname.match(
          /\\/time-parallelization\\/chapter-(\\d+)-/,
        )
        return match ? { item, chapter: Number(match[1]) } : null
      })
      .filter((entry) => entry !== null)

    const orderedItems = [...chapterItems].sort((a, b) => a.chapter - b.chapter)
    const needsReordering = orderedItems.some(
      (entry, index) => entry.item !== chapterItems[index]?.item,
    )
    if (needsReordering) {
      for (const { item } of orderedItems) list.appendChild(item)
    }
  }

  document
    .querySelectorAll(".section-ul")
    .forEach((list) => orderTimeParallelChapters(list, ":scope > .section .desc a"))

  const filterExplorer = (explorerRoot) => {
    if (!explorerRoot) return

    for (const list of explorerRoot.querySelectorAll("ul.content")) {
      orderTimeParallelChapters(list, ":scope > a")
    }

    for (const item of explorerRoot.children) {
      if (!(item instanceof HTMLElement) || item.classList.contains("overflow-end")) continue
      const link = item.querySelector(":scope > .tree-item-self a, :scope > a")
      if (!link) continue
      const path = new URL(link.href, window.location.href).pathname
      const englishTarget =
        path === "/en/" ||
        path.startsWith("/en/") ||
        path === basePath + "/en/" ||
        path.startsWith(basePath + "/en/")
      item.toggleAttribute("hidden", isEnglish ? !englishTarget : englishTarget)
    }
  }

  const explorerRoot = document.querySelector(".explorer > .explorer-content > .explorer-ul")
  filterExplorer(explorerRoot)

  if (explorerRoot && !explorerRoot.dataset.languageObserver) {
    const observer = new MutationObserver(() => filterExplorer(explorerRoot))
    observer.observe(explorerRoot, { childList: true, subtree: true })
    explorerRoot.dataset.languageObserver = "true"
    window.addCleanup?.(() => observer.disconnect())
  }

  if (isEnglish) {
    document.querySelectorAll("p").forEach((paragraph) => {
      const match = paragraph.textContent?.trim().match(/^此文件夹下有(\\d+)条笔记。$/)
      if (match) {
        const count = Number(match[1])
        paragraph.textContent = count === 1 ? "1 note in this folder." : count + " notes in this folder."
      }
    })

    document.querySelectorAll("time[datetime]").forEach((time) => {
      const value = time.getAttribute("datetime")
      if (!value) return
      const date = new Date(value)
      if (!Number.isNaN(date.valueOf())) {
        time.textContent = date.toLocaleDateString("en-US", {
          year: "numeric",
          month: "short",
          day: "numeric",
        })
      }
    })

    const contentMeta = document.querySelector(".content-meta")
    if (contentMeta) {
      for (const span of contentMeta.querySelectorAll("span")) {
        if (span.textContent?.includes("分钟阅读")) {
          span.textContent = span.textContent.replace(/(\\d+)分钟阅读/, "$1 min read")
        }
      }
    }
  }

  if (searchLayout && !searchLayout.dataset.languageObserver) {
    const observer = new MutationObserver(() => filterLanguageLinks(searchLayout))
    observer.observe(searchLayout, { childList: true, subtree: true })
    searchLayout.dataset.languageObserver = "true"
    window.addCleanup?.(() => observer.disconnect())
  }
}

document.addEventListener("nav", localizeKnowledgeBaseUi)
document.addEventListener("render", localizeKnowledgeBaseUi)
localizeKnowledgeBaseUi()
`

export default (() => LanguageSwitcher) satisfies QuartzComponentConstructor
