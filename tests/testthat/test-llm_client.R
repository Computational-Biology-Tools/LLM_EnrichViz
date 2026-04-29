test_that("get_default_model returns correct defaults", {
  expect_equal(get_default_model("gemini"),   "gemini-1.5-flash-latest")
  expect_equal(get_default_model("openai"),   "gpt-4o")
  expect_equal(get_default_model("claude"),   "claude-sonnet-4-6")
  expect_equal(get_default_model("deepseek"), "deepseek-chat")
})

test_that("get_default_model errors on unknown provider", {
  expect_error(get_default_model("unknown_llm"), "Unknown provider")
})

test_that("get_llm_function returns callable functions", {
  for (prov in c("gemini", "openai", "claude", "deepseek")) {
    fn <- get_llm_function(prov)
    expect_true(is.function(fn))
  }
})

test_that("get_llm_function errors on unknown provider", {
  expect_error(get_llm_function("bard"), "Unsupported provider")
})
